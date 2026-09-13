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

PrsmPtrVec readPrsmsWithSpectra(const PrsmParaPtr& prsm_para_ptr,
                                const std::string& input_file_ext) {
  std::string sp_file_name = prsm_para_ptr->getSpectrumFileName();
  std::string input_file_name =
      file_util::basename(sp_file_name) + "." + input_file_ext;
  std::string db_file_name = prsm_para_ptr->getSearchDbFileNameWithFolder();
  FastaIndexReaderPtr seq_reader =
      std::make_shared<FastaIndexReader>(db_file_name);
  ModPtrVec fix_mod_ptr_vec = prsm_para_ptr->getFixModPtrVec();
  PrsmReader prsm_reader(input_file_name);
  PrsmPtr prsm_ptr = prsm_reader.readOnePrsm(seq_reader, fix_mod_ptr_vec);

  int group_spec_num = prsm_para_ptr->getGroupSpecNum();
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  MsAlignReaderPtr ms_reader_ptr = std::make_shared<MsAlignReader>(
      sp_file_name, group_spec_num, sp_para_ptr->getActivationPtr());
  DeconvMsPtrVec deconv_ms_ptr_vec = ms_reader_ptr->getNextMsPtrVec();
  PrsmPtrVec prsm_list;
  while (deconv_ms_ptr_vec.size() != 0) {
    MsHeaderPtr header_ptr = deconv_ms_ptr_vec[0]->getMsHeaderPtr();
    if (header_ptr->containsPrec()) {
      double prec_mono_mass =
          header_ptr->getFirstPrecMonoMass() - sp_para_ptr->getNTermLabelMass();
      SpectrumSetPtr spec_set_ptr = spectrum_set_factory::geneSpectrumSetPtr(
          deconv_ms_ptr_vec, sp_para_ptr, prec_mono_mass);
      if (spec_set_ptr->isValid()) {
        int spec_id = spec_set_ptr->getSpectrumId();
        while (prsm_ptr != nullptr && prsm_ptr->getSpectrumId() == spec_id) {
          DeconvMsPtrVec set_ms_ptr_vec = spec_set_ptr->getDeconvMsPtrVec();
          prsm_ptr->setDeconvMsPtrVec(set_ms_ptr_vec);
          double new_prec_mass = prsm_ptr->getAdjustedPrecMass();
          ExtendMsPtrVec extend_ms_ptr_vec =
              extend_ms_factory::geneMsThreePtrVec(set_ms_ptr_vec, sp_para_ptr,
                                                   new_prec_mass);
          prsm_ptr->setRefineMsVec(extend_ms_ptr_vec);
          prsm_list.push_back(prsm_ptr);
          prsm_ptr = prsm_reader.readOnePrsm(seq_reader, fix_mod_ptr_vec);
        }
      }
    }
    deconv_ms_ptr_vec = ms_reader_ptr->getNextMsPtrVec();
  }
  prsm_reader.close();
  return prsm_list;
}

}  // namespace prsm_reader_util

}  // namespace toppic
