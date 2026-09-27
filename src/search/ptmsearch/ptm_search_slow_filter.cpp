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

#include "search/ptmsearch/ptm_search_slow_filter.hpp"

#include <string>

#include "common/base/prot_mod_base.hpp"
#include "common/util/file_util.hpp"
#include "seq/proteoform_factory.hpp"

namespace toppic {

PtmSearchSlowFilter::PtmSearchSlowFilter(
    const SpectrumSetPtr& spectrum_set_ptr,
    const SimplePrsmPtrVec& simple_prsm_ptrs, const PtmSearchMngPtr& mng_ptr) {
  std::string db_file_name =
      mng_ptr->prsm_para_ptr_->getSearchDbFileNameWithFolder();
  // one indexed FASTA reader per worker thread: a reader was constructed -
  // an fai_load parsing the 40k-entry .fai index - for every spectrum task
  static thread_local FastaIndexReaderPtr reader_ptr;
  if (reader_ptr == nullptr) {
    reader_ptr = std::make_shared<FastaIndexReader>(db_file_name);
  }
  ModPtrVec fix_mod_list = mng_ptr->prsm_para_ptr_->getFixModPtrVec();
  ProtModPtrVec prot_mod_ptr_vec = mng_ptr->prsm_para_ptr_->getProtModPtrVec();
  ProteoformPtrVec proteoform_ptr_vec;
  for (size_t i = 0; i < simple_prsm_ptrs.size(); i++) {
    std::string seq_name = simple_prsm_ptrs[i]->getSeqName();
    std::string seq_desc = simple_prsm_ptrs[i]->getSeqDesc();
    ProteoformPtr proteo_ptr = proteoform_factory::readFastaToProteoformPtr(
        reader_ptr, seq_name, seq_desc, fix_mod_list);
    ProteoformPtrVec mod_form_ptr_vec =
        proteoform_factory::geneProtModProteoform(proteo_ptr, prot_mod_ptr_vec);
    for (size_t j = 0; j < mod_form_ptr_vec.size(); j++) {
      // this slow match is used for four type of alignment, INTERNAL is used to
      // include more diagonals.
      PtmSlowMatchPtr ptm_slow_match_ptr =
          std::make_shared<PtmSlowMatch>(mod_form_ptr_vec[j], spectrum_set_ptr,
                                         ProteoformType::INTERNAL, mng_ptr);
      slow_match_ptrs_.push_back(ptm_slow_match_ptr);
    }
  }

  // Compute complete, prefix, suffix and internal PrSMs together, one fused
  // dynamic program per candidate (search/ptmsearch/ps_align_multi_type.hpp)
  // instead of up to four separate PsAlign passes. Suffix and internal
  // widen the search to fragments that do not start/end at either end of
  // the full-length protein, so, as before, they are only computed for
  // proteoforms without an N/C-terminal protein modification.
  for (size_t i = 0; i < slow_match_ptrs_.size(); i++) {
    bool is_none_mod = slow_match_ptrs_[i]->getProteoform()->getProtModPtr() ==
                       ProtModBase::getProtModPtr_NONE();
    ProteoformTypePtrVec types;
    if (mng_ptr->prsm_para_ptr_->allowCompleteProt()) {
      types.push_back(ProteoformType::COMPLETE);
    }
    if (mng_ptr->prsm_para_ptr_->allowPrefixProt()) {
      types.push_back(ProteoformType::PREFIX);
    }
    if (is_none_mod && mng_ptr->prsm_para_ptr_->allowSuffixProt()) {
      types.push_back(ProteoformType::SUFFIX);
    }
    if (is_none_mod && mng_ptr->prsm_para_ptr_->allowCompleteProt()) {
      // Preserves the original condition (also gated on allowCompleteProt,
      // not allowInternalProt) so results match the existing implementation.
      types.push_back(ProteoformType::INTERNAL);
    }

    std::vector<PrsmPtrVec> result_2d;
    if (!types.empty()) {
      slow_match_ptrs_[i]->compute(types, result_2d);
    }

    size_t idx = 0;
    complete_prsm_2d_ptrs_.push_back(
        mng_ptr->prsm_para_ptr_->allowCompleteProt() ? result_2d[idx++]
                                                     : PrsmPtrVec());
    prefix_prsm_2d_ptrs_.push_back(mng_ptr->prsm_para_ptr_->allowPrefixProt()
                                       ? result_2d[idx++]
                                       : PrsmPtrVec());
    if (is_none_mod) {
      suffix_prsm_2d_ptrs_.push_back(mng_ptr->prsm_para_ptr_->allowSuffixProt()
                                         ? result_2d[idx++]
                                         : PrsmPtrVec());
      internal_prsm_2d_ptrs_.push_back(
          mng_ptr->prsm_para_ptr_->allowCompleteProt() ? result_2d[idx++]
                                                       : PrsmPtrVec());
    }
  }
}

PrsmPtrVec PtmSearchSlowFilter::getPrsms(int shift_num,
                                         const PrsmParaPtr& prsm_para_ptr,
                                         const ProteoformTypePtr& type_ptr) {
  PrsmPtrVec prsm_ptrs;
  if (type_ptr == ProteoformType::COMPLETE &&
      prsm_para_ptr->allowCompleteProt()) {
    for (size_t i = 0; i < complete_prsm_2d_ptrs_.size(); i++) {
      if (complete_prsm_2d_ptrs_[i][shift_num] != nullptr) {
        prsm_ptrs.push_back(complete_prsm_2d_ptrs_[i][shift_num]);
      }
    }
  } else if (type_ptr == ProteoformType::PREFIX &&
             prsm_para_ptr->allowPrefixProt()) {
    for (size_t i = 0; i < prefix_prsm_2d_ptrs_.size(); i++) {
      if (prefix_prsm_2d_ptrs_[i][shift_num] != nullptr) {
        prsm_ptrs.push_back(prefix_prsm_2d_ptrs_[i][shift_num]);
      }
    }
  } else if (type_ptr == ProteoformType::SUFFIX &&
             prsm_para_ptr->allowSuffixProt()) {
    for (size_t i = 0; i < suffix_prsm_2d_ptrs_.size(); i++) {
      if (suffix_prsm_2d_ptrs_[i][shift_num] != nullptr) {
        prsm_ptrs.push_back(suffix_prsm_2d_ptrs_[i][shift_num]);
      }
    }
  } else if (type_ptr == ProteoformType::INTERNAL &&
             prsm_para_ptr->allowInternalProt()) {
    for (size_t i = 0; i < internal_prsm_2d_ptrs_.size(); i++) {
      if (internal_prsm_2d_ptrs_[i][shift_num] != nullptr) {
        prsm_ptrs.push_back(internal_prsm_2d_ptrs_[i][shift_num]);
      }
    }
  }
  return prsm_ptrs;
}

} /* namespace toppic */
