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

#include "topfd/deconv/deconv_ms3_process.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "common/util/logger.hpp"
#include "ms/env/avg_mass_util.hpp"
#include "ms/env/match_env_util.hpp"
#include "ms/mzml/pw_ms_reader.hpp"
#include "ms/spec/msalign_reader.hpp"
#include "ms/spec/msalign_writer.hpp"
#include "ms/spec/peak_util.hpp"
#include "ms/spec/precursor.hpp"
#include "topfd/deconv/deconv_single_sp.hpp"
#include "topfd/deconv/greedy_ms2_deconv.hpp"

namespace toppic {

namespace {

// The isotope peak spacing in m/z at charge z
constexpr double ISOTOPE_DIFF = 1.00235;

// The MS2 species whose isotopic envelope has a peak in the MS3 isolation
// window [win_begin, win_end]: the most intense one, or -1. An envelope's
// peaks are the monoisotopic m/z and the 1 + mass / 800 isotope peaks above
// it (the heavier the species, the more peaks carry intensity).
int findPrecSpecies(const DeconvMsPtr& ms2_ptr, double win_begin,
                    double win_end) {
  int best = -1;
  double best_inte = -1;
  for (size_t i = 0; i < ms2_ptr->size(); i++) {
    DeconvPeakPtr peak_ptr = ms2_ptr->getPeakPtr(i);
    double mass = peak_ptr->getPosition();
    int charge = peak_ptr->getCharge();
    if (charge < 1) continue;
    double mono_mz = peak_util::compMz(mass, charge);
    int peak_num = 2 + static_cast<int>(mass / 800);
    bool in_window = false;
    for (int k = 0; k < peak_num; k++) {
      double mz = mono_mz + k * ISOTOPE_DIFF / charge;
      if (mz >= win_begin && mz <= win_end) {
        in_window = true;
        break;
      }
    }
    if (in_window && peak_ptr->getIntensity() > best_inte) {
      best = static_cast<int>(i);
      best_inte = peak_ptr->getIntensity();
    }
  }
  return best;
}

}  // namespace

DeconvMs3Process::DeconvMs3Process(const TopfdParaPtr& topfd_para_ptr)
    : topfd_para_ptr_(topfd_para_ptr) {}

void DeconvMs3Process::readMsAlign(const std::string& file_name,
                                   std::map<int, DeconvMsPtr>& ms_map) {
  if (!std::ifstream(file_name).good()) {
    LOG_WARN("The msalign file " << file_name << " does not exist.");
    return;
  }
  MsAlignReader reader(file_name);
  DeconvMsPtr ms_ptr = reader.getNextMsPtr();
  while (ms_ptr != nullptr) {
    ms_map[ms_ptr->getMsHeaderPtr()->getFirstScanNum()] = ms_ptr;
    ms_ptr = reader.getNextMsPtr();
  }
}

void DeconvMs3Process::process() {
  std::string output_base_name = topfd_para_ptr_->getOutputBaseName();
  std::map<int, DeconvMsPtr> ms1_map;
  std::map<int, DeconvMsPtr> ms2_map;
  readMsAlign(output_base_name + "_ms1.msalign", ms1_map);
  readMsAlign(output_base_name + "_ms2.msalign", ms2_map);

  PwMsReader reader(topfd_para_ptr_->getMzmlFileName(),
                    topfd_para_ptr_->getPrecWindowWidth(),
                    topfd_para_ptr_->getActivation());
  std::string msalign_name = output_base_name + "_ms3.msalign";
  MsAlignWriter writer(msalign_name);
  writer.writePara(topfd_para_ptr_->getParaStr("#", "\t"));
  std::ofstream tsv(output_base_name + "_ms3.tsv");
  tsv << "Spectrum ID\tScan\tMS2 scan\tMS2 ID\tMS1 scan\tMS1 ID\tRetention time"
      << "\tMS2 precursor m/z\tMS2 precursor charge\tMS2 precursor mass"
      << "\tIsolation target m/z\tIsolation window begin\tIsolation window end"
      << "\tPrecursor source\tPrecursor m/z\tPrecursor charge\tPrecursor mass"
      << "\tPrecursor intensity\tPeaks\tDeconvoluted masses" << std::endl;

  int total = topfd_para_ptr_->getMs3ScanNum();
  int spec_id = 0;
  int species_num = 0;
  int selected_num = 0;
  int no_parent_num = 0;
  int cnt = 0;
  while (reader.readNext() != -1) {
    MsHeaderPtr header_ptr = reader.getHeaderPtr();
    if (header_ptr == nullptr || header_ptr->getMsLevel() != 3 ||
        header_ptr->isPtcr() || header_ptr->isFragPtcr() ||
        header_ptr->isReducedPrec()) {
      continue;
    }
    cnt++;
    if (cnt % 100 == 0 || cnt == total) {
      std::cout << "\r" << std::flush << "Processing MS3 spectrum scan "
                << header_ptr->getFirstScanNum() << " ... "
                << (total > 0 ? cnt * 100 / total : 0) << "% finished.";
    }
    PeakPtrVec peak_list = reader.getPeakList();
    // the parent MS2 spectrum and the cycle's MS1 spectrum
    int ms2_scan = header_ptr->getParentScan();
    auto ms2_it = ms2_map.find(ms2_scan);
    DeconvMsPtr ms2_ptr = ms2_it == ms2_map.end() ? nullptr : ms2_it->second;
    auto ms1_it = ms1_map.find(header_ptr->getMsOneScan());
    header_ptr->setMsOneId(ms1_it == ms1_map.end()
                               ? -1
                               : ms1_it->second->getMsHeaderPtr()->getSpecId());
    double ms2_prec_mz = -1;
    int ms2_prec_charge = 0;
    double ms2_prec_mass = -1;
    if (ms2_ptr != nullptr) {
      header_ptr->setMsTwoId(ms2_ptr->getMsHeaderPtr()->getSpecId());
      const PrecursorPtrVec& ms2_precs =
          ms2_ptr->getMsHeaderPtr()->getPrecPtrVec();
      if (!ms2_precs.empty()) {
        ms2_prec_mz = ms2_precs[0]->getMonoMz();
        ms2_prec_charge = ms2_precs[0]->getCharge();
        ms2_prec_mass = ms2_precs[0]->getMonoMass();
      }
    } else {
      no_parent_num++;
    }
    header_ptr->setMsTwoPrec(ms2_prec_mz, ms2_prec_charge, ms2_prec_mass);
    // the precursor: the MS2 species in the isolation window, else the
    // selected ion at every charge up to the MS2 precursor's
    double win_begin = header_ptr->getPrecWinBegin();
    double win_end = header_ptr->getPrecWinEnd();
    PrecursorPtrVec prec_ptr_vec;
    std::string source;
    double max_mass = topfd_para_ptr_->getMaxMass();
    int max_charge = topfd_para_ptr_->getMaxCharge();
    int species = ms2_ptr == nullptr
                      ? -1
                      : findPrecSpecies(ms2_ptr, win_begin, win_end);
    if (species >= 0) {
      DeconvPeakPtr peak_ptr = ms2_ptr->getPeakPtr(species);
      double mass = peak_ptr->getPosition();
      int charge = peak_ptr->getCharge();
      prec_ptr_vec.push_back(std::make_shared<Precursor>(
          0, -1, peak_util::compMz(mass, charge), charge,
          peak_ptr->getIntensity()));
      max_mass = std::min(max_mass, mass + 1.0);
      max_charge = std::min(max_charge, charge);
      source = "ms2_species";
      species_num++;
    } else {
      double sel_mz = header_ptr->getSelectedMz() > 0
                          ? header_ptr->getSelectedMz()
                          : header_ptr->getPrecTargetMz();
      int top_charge = ms2_prec_charge > 0 ? ms2_prec_charge : max_charge;
      top_charge = std::min(top_charge, max_charge);
      for (int z = 1; z <= top_charge; z++) {
        prec_ptr_vec.push_back(
            std::make_shared<Precursor>(z - 1, -1, sel_mz, z, 0.0));
      }
      if (ms2_prec_mass > 0) max_mass = std::min(max_mass, ms2_prec_mass + 1.0);
      max_charge = top_charge;
      source = "selected_ion";
      selected_num++;
    }
    // the deconvolution, the MS/MS method bounded by the precursor
    MatchEnvPtrVec deconv_envs;
    if (!peak_list.empty() && topfd_para_ptr_->isMs2GreedyDeconv()) {
      deconv_envs = greedy_ms2_deconv::deconv(peak_list, topfd_para_ptr_,
                                              max_mass, max_charge, "");
    } else if (!peak_list.empty()) {
      DeconvSingleSpPtr deconv_ptr = std::make_shared<DeconvSingleSp>(
          topfd_para_ptr_, peak_list, 2, max_mass, max_charge);
      deconv_envs = deconv_ptr->deconv();
    }
    header_ptr->setSpecId(spec_id++);
    header_ptr->setPrecPtrVec(prec_ptr_vec);
    DeconvMsPtr deconv_ms_ptr =
        match_env_util::getDeconvMsPtr(header_ptr, deconv_envs);
    if (topfd_para_ptr_->isOutputAvgMass()) {
      deconv_ms_ptr = avg_mass_util::getAvgMassMs(deconv_ms_ptr);
    }
    writer.writeMs(deconv_ms_ptr);
    tsv << std::fixed << header_ptr->getSpecId() << "\t"
        << header_ptr->getFirstScanNum() << "\t" << ms2_scan << "\t"
        << header_ptr->getMsTwoId() << "\t" << header_ptr->getMsOneScan()
        << "\t" << header_ptr->getMsOneId() << "\t" << std::setprecision(3)
        << header_ptr->getRetentionTime() / 60 << "\t" << std::setprecision(5)
        << ms2_prec_mz << "\t" << ms2_prec_charge << "\t" << ms2_prec_mass
        << "\t" << header_ptr->getPrecTargetMz() << "\t" << win_begin << "\t"
        << win_end << "\t" << source << "\t";
    for (size_t i = 0; i < prec_ptr_vec.size(); i++) {
      tsv << (i > 0 ? ":" : "") << prec_ptr_vec[i]->getMonoMz();
    }
    tsv << "\t";
    for (size_t i = 0; i < prec_ptr_vec.size(); i++) {
      tsv << (i > 0 ? ":" : "") << prec_ptr_vec[i]->getCharge();
    }
    tsv << "\t";
    for (size_t i = 0; i < prec_ptr_vec.size(); i++) {
      tsv << (i > 0 ? ":" : "") << prec_ptr_vec[i]->getMonoMass();
    }
    tsv << "\t" << std::setprecision(2)
        << (prec_ptr_vec.empty() ? 0.0 : prec_ptr_vec[0]->getInte()) << "\t"
        << peak_list.size() << "\t" << deconv_ms_ptr->size() << std::endl;
  }
  std::cout << std::endl;
  std::cout << "MS3 spectra: " << cnt << " deconvoluted, " << species_num
            << " with the precursor from an MS2 species, " << selected_num
            << " from the selected ion, " << no_parent_num
            << " without a parent MS2 spectrum." << std::endl;
}

}  // namespace toppic
