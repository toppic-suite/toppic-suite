// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "prsm/prsm_second_prec.hpp"

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <vector>

#include "common/base/mass_constant.hpp"
#include "ms/spec/deconv_ms.hpp"
#include "ms/spec/msalign_reader.hpp"
#include "ms/spec/msalign_util.hpp"
#include "ms/spec/msalign_writer.hpp"
#include "ms/spec/precursor.hpp"
#include "para/sp_para.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_reader_util.hpp"
#include "prsm/prsm_str.hpp"
#include "prsm/prsm_xml_writer.hpp"
#include "prsm/theo_peak_util.hpp"
#include "seq/fasta_index_reader.hpp"

namespace toppic {

namespace prsm_second_prec {

namespace {

// the ids of the spectrum's masses within the fragment tolerance of one of
// the proteoform's fragment masses, as it is or one isotope lower or higher
std::set<int> matchedPeakIds(const ProteoformPtr& proteoform_ptr,
                             const DeconvMsPtr& ms_ptr,
                             const SpParaPtr& sp_para_ptr) {
  MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
  TheoPeakPtrVec theo_peaks = theo_peak_util::geneProteoformTheoPeak(
      proteoform_ptr, header_ptr->getActivationPtr(),
      sp_para_ptr->getMinMass(), header_ptr->getNTermLabelMass());
  std::vector<double> theo_masses = theo_peak_util::getTheoMassVec(theo_peaks);
  std::sort(theo_masses.begin(), theo_masses.end());
  PeakTolerancePtr tole_ptr = sp_para_ptr->getFragTolerancePtr();
  const double offsets[] = {0.0, -mass_constant::getIsotopeMass(),
                            mass_constant::getIsotopeMass()};
  std::set<int> ids;
  for (const DeconvPeakPtr& peak_ptr : ms_ptr->getPeakPtrVec()) {
    double mass = peak_ptr->getMonoMass();
    double err = tole_ptr->compStrictErrorTole(mass);
    for (double offset : offsets) {
      double m = mass + offset;
      auto it = std::lower_bound(theo_masses.begin(), theo_masses.end(),
                                 m - err);
      if (it != theo_masses.end() && *it <= m + err) {
        ids.insert(peak_ptr->getPeakId());
        break;
      }
    }
  }
  return ids;
}

std::string stripDecoy(const std::string& name) {
  return name.find("DECOY_") == 0 ? name.substr(6) : name;
}

}  // namespace

int writeResidualSpectra(const PrsmParaPtr& prsm_para_ptr,
                         const std::string& record_file_name,
                         const std::string& ori_sp_file_name,
                         const std::string& cutoff_type, double cutoff_value,
                         const std::string& resid_sp_file_name) {
  // the accepted PrSMs of the first pass, by spectrum id (one per spectrum:
  // the first pass searches one precursor)
  std::string db_file_name = prsm_para_ptr->getSearchDbFileNameWithFolder();
  FastaIndexReaderPtr seq_reader =
      std::make_shared<FastaIndexReader>(db_file_name);
  ModPtrVec fix_mod_ptr_vec = prsm_para_ptr->getFixModPtrVec();
  PrsmReader prsm_reader(record_file_name);
  std::map<int, PrsmPtr> accepted;
  PrsmPtr prsm_ptr = prsm_reader.readOnePrsm(seq_reader, fix_mod_ptr_vec);
  bool fdr_cutoff = cutoff_type == "FDR";
  while (prsm_ptr != nullptr) {
    bool decoy =
        prsm_ptr->getProteoformPtr()->getSeqName().find("DECOY_") == 0;
    bool pass = fdr_cutoff ? prsm_ptr->getFdr() >= 0 &&
                                 prsm_ptr->getFdr() <= cutoff_value
                           : prsm_ptr->getEValue() <= cutoff_value;
    if (!decoy && pass) {
      auto it = accepted.find(prsm_ptr->getSpectrumId());
      if (it == accepted.end() ||
          prsm_ptr->getEValue() < it->second->getEValue()) {
        accepted[prsm_ptr->getSpectrumId()] = prsm_ptr;
      }
    }
    prsm_ptr = prsm_reader.readOnePrsm(seq_reader, fix_mod_ptr_vec);
  }
  prsm_reader.close();

  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  // every mass of the spectra, as the post mass matching reads them
  MsAlignReader reader(ori_sp_file_name, 1, sp_para_ptr->getActivationPtr());
  MsAlignWriter writer(resid_sp_file_name);
  writer.writePara(msalign_util::readParaStr(ori_sp_file_name));
  int written = 0;
  int removed_spec_num = 0;
  int removed_mass_num = 0;
  DeconvMsPtr ms_ptr = reader.getNextMsPtr();
  while (ms_ptr != nullptr) {
    MsHeaderPtr header_ptr = ms_ptr->getMsHeaderPtr();
    const PrecursorPtrVec& precs = header_ptr->getPrecPtrVec();
    if (precs.size() >= 2) {
      // the second precursor first, then the first and the others, with
      // new precursor ids
      MsHeaderPtr new_header_ptr = std::make_shared<MsHeader>(*header_ptr);
      std::vector<size_t> order = {1, 0};
      for (size_t i = 2; i < precs.size(); i++) {
        order.push_back(i);
      }
      PrecursorPtrVec new_precs;
      for (size_t k = 0; k < order.size(); k++) {
        const PrecursorPtr& p = precs[order[k]];
        new_precs.push_back(std::make_shared<Precursor>(
            static_cast<int>(k), p->getFeatureId(), p->getMonoMz(),
            p->getCharge(), p->getInte()));
      }
      new_header_ptr->setPrecPtrVec(new_precs);
      DeconvPeakPtrVec peaks = ms_ptr->getPeakPtrVec();
      auto it = accepted.find(header_ptr->getSpecId());
      if (it != accepted.end()) {
        std::set<int> ids = matchedPeakIds(it->second->getProteoformPtr(),
                                           ms_ptr, sp_para_ptr);
        if (!ids.empty()) {
          DeconvPeakPtrVec kept;
          for (const DeconvPeakPtr& peak_ptr : peaks) {
            if (ids.count(peak_ptr->getPeakId()) == 0) {
              kept.push_back(peak_ptr);
            }
          }
          removed_spec_num++;
          removed_mass_num += static_cast<int>(peaks.size() - kept.size());
          peaks = kept;
        }
      }
      writer.writeMs(std::make_shared<DeconvMs>(new_header_ptr, peaks));
      written++;
    }
    ms_ptr = reader.getNextMsPtr();
  }
  std::cout << "Second precursor search - " << written
            << " spectra with a second precursor, " << removed_mass_num
            << " masses of " << removed_spec_num
            << " identified spectra removed." << std::endl;
  return written;
}

void mergeSpectra(const std::string& first_sp_file_name,
                  const std::string& second_sp_file_name,
                  const std::string& out_sp_file_name) {
  std::ofstream outfile(out_sp_file_name);
  outfile << msalign_util::readParaStr(first_sp_file_name) << std::endl;
  const std::string names[] = {first_sp_file_name, second_sp_file_name};
  for (int i = 0; i < 2; i++) {
    int offset = i * SpPara::getMaxSpecNumPerFile();
    MsAlignReader sp_reader(names[i]);
    std::vector<std::string> ms_lines = sp_reader.readOneStrSpectrum();
    while (ms_lines.size() > 0) {
      for (size_t k = 0; k < ms_lines.size(); k++) {
        if (ms_lines[k].substr(0, 12) == "SPECTRUM_ID=") {
          outfile << "SPECTRUM_ID="
                  << (offset + std::stoi(ms_lines[k].substr(12)))
                  << std::endl;
        } else if (ms_lines[k].substr(0, 10) == "MS_ONE_ID=") {
          outfile << "MS_ONE_ID="
                  << (offset + std::stoi(ms_lines[k].substr(10)))
                  << std::endl;
        } else {
          outfile << ms_lines[k] << std::endl;
        }
      }
      outfile << std::endl;
      ms_lines = sp_reader.readOneStrSpectrum();
    }
  }
  outfile.close();
}

void mergePrsms(const std::string& first_prsm_file_name,
                const std::string& second_prsm_file_name,
                const std::string& out_prsm_file_name,
                const std::string& sp_short_name) {
  PrsmStrPtrVec first = prsm_reader_util::readAllPrsmStrs(first_prsm_file_name);
  PrsmStrPtrVec second =
      prsm_reader_util::readAllPrsmStrs(second_prsm_file_name);
  // the units of every spectrum, the first pass's first
  std::map<int, PrsmStrPtrVec> by_spec;
  for (const PrsmStrPtr& p : first) {
    by_spec[p->getSpectrumId()].push_back(p);
  }
  int kept = 0;
  int replaced = 0;
  int dropped = 0;
  int offset = SpPara::getMaxSpecNumPerFile();
  for (const PrsmStrPtr& p : second) {
    p->setSpectrumId(p->getSpectrumId() + offset);
    p->setPrecursorId(1);
    p->setFileName(sp_short_name);
    PrsmStrPtrVec& units = by_spec[p->getSpectrumId() - offset];
    bool same = false;
    for (size_t i = 0; i < units.size(); i++) {
      const PrsmStrPtr& q = units[i];
      if (stripDecoy(q->getSeqName()) == stripDecoy(p->getSeqName()) ||
          (q->getFracFeatureId() >= 0 &&
           q->getFracFeatureId() == p->getFracFeatureId())) {
        same = true;
        if (p->getEValue() < q->getEValue()) {
          units[i] = p;
          replaced++;
        } else {
          dropped++;
        }
        break;
      }
    }
    if (!same) {
      units.push_back(p);
      kept++;
    }
  }
  PrsmStrPtrVec all;
  for (auto& spec : by_spec) {
    all.insert(all.end(), spec.second.begin(), spec.second.end());
  }
  std::sort(all.begin(), all.end(), PrsmStr::cmpSpecIncPrecIncEvalueIncProtInc);
  PrsmXmlWriter writer(out_prsm_file_name);
  writer.writeVector(all);
  writer.close();
  std::cout << "Second precursor search - " << first.size()
            << " first-pass PrSMs; of " << second.size()
            << " second-pass PrSMs " << kept << " added, " << replaced
            << " replacing a first-pass PrSM of the same protein or feature, "
            << dropped << " dropped for one." << std::endl;
}

}  // namespace prsm_second_prec

}  // namespace toppic
