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

#include "prsm/prsm_reader_util.hpp"

#include <cstddef>
#include <map>
#include <vector>

#include "common/util/file_util.hpp"
#include "common/util/logger.hpp"
#include "ms/factory/extend_ms_factory.hpp"
#include "ms/factory/spectrum_set_factory.hpp"
#include "ms/spec/msalign_reader.hpp"
#include "prsm/prsm_reader.hpp"

namespace toppic {

namespace prsm_reader_util {

PrsmStrPtrVec readAllPrsmStrs(const std::string& input_file_name) {
  PrsmReader reader(input_file_name);
  PrsmStrPtrVec prsm_str_ptrs;
  PrsmStrPtr prsm_str_ptr = reader.readOnePrsmStr();
  while (prsm_str_ptr != nullptr) {
    prsm_str_ptrs.push_back(prsm_str_ptr);
    prsm_str_ptr = reader.readOnePrsmStr();
  }
  reader.close();
  return prsm_str_ptrs;
}

PrsmStrPtrVec readAllPrsmStrsMatchSeq(const std::string& input_file_name) {
  PrsmReaderPtr str_reader = std::make_shared<PrsmReader>(input_file_name);
  PrsmStrPtrVec prsm_str_ptrs;
  PrsmStrPtr prsm_str_ptr = str_reader->readOnePrsmStr();
  while (prsm_str_ptr != nullptr) {
    prsm_str_ptrs.push_back(prsm_str_ptr);
    prsm_str_ptr = str_reader->readOnePrsmStr();
  }
  str_reader->close();
  return prsm_str_ptrs;
}

PrsmPtrVec readAllPrsms(const std::string& prsm_file_name,
                        const std::string& db_file_name,
                        const ModPtrVec& fix_mod_list) {
  FastaIndexReaderPtr fasta_reader_ptr =
      std::make_shared<FastaIndexReader>(db_file_name);
  PrsmReader reader(prsm_file_name);
  PrsmPtrVec prsm_ptrs;
  PrsmPtr prsm_ptr = reader.readOnePrsm(fasta_reader_ptr, fix_mod_list);
  while (prsm_ptr != nullptr) {
    prsm_ptrs.push_back(prsm_ptr);
    prsm_ptr = reader.readOnePrsm(fasta_reader_ptr, fix_mod_list);
  }
  reader.close();
  return prsm_ptrs;
}

PrsmPtrVec readAllPrsms(const std::string& prsm_file_name,
                        const FastaIndexReaderPtr& fasta_reader_ptr,
                        const ModPtrVec& fix_mod_list) {
  PrsmReader reader(prsm_file_name);
  PrsmPtrVec prsm_ptrs;
  PrsmPtr prsm_ptr = reader.readOnePrsm(fasta_reader_ptr, fix_mod_list);
  while (prsm_ptr != nullptr) {
    prsm_ptrs.push_back(prsm_ptr);
    prsm_ptr = reader.readOnePrsm(fasta_reader_ptr, fix_mod_list);
  }
  reader.close();
  return prsm_ptrs;
}

namespace {

// Attaches the spectra of sp_file_name to the PrSMs of their spectrum id
// (by_spec: spectrum id -> the PrSMs' indices in prsm_ptrs), every PrSM of
// the id or, with match_prec, only those whose precursor id is the
// spectrum's first precursor's; the spectrum set is the first precursor's,
// the refined spectra the PrSM's adjusted precursor mass's. attached marks
// the PrSMs that got spectra.
void attachSpectra(const PrsmPtrVec& prsm_ptrs,
                   const std::map<int, std::vector<size_t>>& by_spec,
                   const std::string& sp_file_name,
                   const PrsmParaPtr& prsm_para_ptr, double env_cnn_cutoff,
                   bool match_prec, std::vector<bool>& attached) {
  int group_spec_num = prsm_para_ptr->getGroupSpecNum();
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  MsAlignReaderPtr ms_reader_ptr = std::make_shared<MsAlignReader>(
      sp_file_name, group_spec_num, sp_para_ptr->getActivationPtr(),
      env_cnn_cutoff);
  DeconvMsPtrVec deconv_ms_ptr_vec = ms_reader_ptr->getNextMsPtrVec();
  while (deconv_ms_ptr_vec.size() != 0) {
    MsHeaderPtr header_ptr = deconv_ms_ptr_vec[0]->getMsHeaderPtr();
    auto it = by_spec.find(header_ptr->getSpecId());
    if (it != by_spec.end() && header_ptr->containsPrec()) {
      // the observed precursor mass: the spectrum set subtracts the N-terminal
      // label (toppic -N) itself, for its PRM and extended spectra
      double prec_mono_mass = header_ptr->getFirstPrecMonoMass();
      SpectrumSetPtr spec_set_ptr = spectrum_set_factory::geneSpectrumSetPtr(
          deconv_ms_ptr_vec, sp_para_ptr, prec_mono_mass);
      if (spec_set_ptr->isValid()) {
        int prec_id = header_ptr->getFirstPrecId();
        for (size_t i : it->second) {
          const PrsmPtr& prsm_ptr = prsm_ptrs[i];
          if (match_prec && prsm_ptr->getPrecursorId() != prec_id) {
            continue;
          }
          DeconvMsPtrVec set_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();
          prsm_ptr->setDeconvMsPtrVec(set_ms_ptr_vec);
          double new_prec_mass = prsm_ptr->getAdjustedPrecMass();
          ExtendMsPtrVec extend_ms_ptr_vec =
              extend_ms_factory::geneMsThreePtrVec(set_ms_ptr_vec, sp_para_ptr,
                                                   new_prec_mass);
          prsm_ptr->setRefineMsVec(extend_ms_ptr_vec);
          attached[i] = true;
        }
      }
    }
    deconv_ms_ptr_vec = ms_reader_ptr->getNextMsPtrVec();
  }
}

}  // namespace

PrsmPtrVec readPrsmsWithSpectra(const PrsmParaPtr& prsm_para_ptr,
                                const std::string& input_file_ext,
                                double env_cnn_cutoff,
                                const std::string& second_sp_file_name) {
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string input_file_name =
      file_util::basename(sp_file_name) + "." + input_file_ext;
  std::string db_file_name = prsm_para_ptr->getSearchDbFileNameWithFolder();
  FastaIndexReaderPtr seq_reader =
      std::make_shared<FastaIndexReader>(db_file_name);
  PrsmPtrVec prsm_ptrs = readAllPrsms(input_file_name, seq_reader,
                                      prsm_para_ptr->getFixModPtrVec());
  std::map<int, std::vector<size_t>> by_spec;
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    by_spec[prsm_ptrs[i]->getSpectrumId()].push_back(i);
  }
  std::vector<bool> attached(prsm_ptrs.size(), false);
  bool two_files = second_sp_file_name != "";
  attachSpectra(prsm_ptrs, by_spec, sp_file_name, prsm_para_ptr,
                env_cnn_cutoff, two_files, attached);
  if (two_files) {
    attachSpectra(prsm_ptrs, by_spec, second_sp_file_name, prsm_para_ptr,
                  env_cnn_cutoff, true, attached);
  }
  PrsmPtrVec prsm_list;
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    if (attached[i]) {
      prsm_list.push_back(prsm_ptrs[i]);
    }
  }
  return prsm_list;
}

}  // namespace prsm_reader_util

}  // namespace toppic
