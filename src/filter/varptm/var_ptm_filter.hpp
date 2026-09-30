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

#ifndef TOPPIC_FILTER_VARPTM_VAR_PTM_FILTER_HPP_
#define TOPPIC_FILTER_VARPTM_VAR_PTM_FILTER_HPP_

#include "filter/massmatch/filter_spec.hpp"
#include "filter/massmatch/mass_match.hpp"
#include "filter/massmatch/mass_match_util.hpp"
#include "filter/mng/var_ptm_filter_mng.hpp"
#include "ms/spec/prm_ms.hpp"
#include "prsm/simple_prsm.hpp"

namespace toppic {

class VarPtmFilter {
 public:
  VarPtmFilter(const ProteoformPtrVec& proteo_ptrs,
               const VarPtmFilterMngPtr& mng_ptr, std::string block_str);

  // the block independent part of a spectrum set, computed once per stage
  static FilterSpecPtr geneFilterSpec(const PrmMsPtrVec& prm_ms_ptr_vec,
                                     const PrmMsPtrVec& srm_ms_ptr_vec,
                                     const VarPtmFilterMngPtr& mng_ptr);

  void computeBestMatch(const FilterSpecPtr& spec_ptr);

  void computeBestMatch(const PrmMsPtrVec& prm_ms_ptr_vec,
                        const PrmMsPtrVec& srm_ms_ptr_vec);

  SimplePrsmPtrVec getCompMatchPtrs() { return comp_match_ptrs_; }
  SimplePrsmPtrVec getPrefMatchPtrs() { return pref_match_ptrs_; }
  SimplePrsmPtrVec getSuffMatchPtrs() { return suff_match_ptrs_; }
  SimplePrsmPtrVec getInternalMatchPtrs() { return internal_match_ptrs_; }

  MassMatchPtr getDiagIndexPtr() { return diag_index_ptr_; }
  MassMatchPtr getRevDiagIndexPtr() { return rev_diag_index_ptr_; }
  MassMatchPtr getTermIndexPtr() { return term_index_ptr_; }
  MassMatchPtr getRevTermIndexPtr() { return rev_term_index_ptr_; }

 private:
  VarPtmFilterMngPtr mng_ptr_;
  ProteoformPtrVec proteo_ptrs_;

  // Score arrays of the four indexes, the rows of each that received the
  // precursor bonus, and those rows grouped by proteoform, which the four
  // proteoform types share. They are members so that they are allocated once
  // per database block instead of once per spectrum.
  std::vector<short> term_scores_;
  std::vector<short> diag_scores_;
  std::vector<short> rev_term_scores_;
  std::vector<short> rev_diag_scores_;
  std::vector<int> term_prec_rows_;
  std::vector<int> diag_prec_rows_;
  std::vector<int> rev_term_prec_rows_;
  std::vector<int> rev_diag_prec_rows_;
  mass_match_util::PrecMatchRows term_rows_2d_;
  mass_match_util::PrecMatchRows diag_rows_2d_;
  mass_match_util::PrecMatchRows rev_term_rows_2d_;
  mass_match_util::PrecMatchRows rev_diag_rows_2d_;

  MassMatchPtr diag_index_ptr_;
  MassMatchPtr rev_diag_index_ptr_;
  MassMatchPtr term_index_ptr_;
  MassMatchPtr rev_term_index_ptr_;

  // the sequence ranges of the index rows
  mass_match_util::IndexRowRanges row_ranges_;

  SimplePrsmPtrVec comp_match_ptrs_;
  SimplePrsmPtrVec pref_match_ptrs_;
  SimplePrsmPtrVec suff_match_ptrs_;
  SimplePrsmPtrVec internal_match_ptrs_;

  // which db_block it is. Used when merging the output files by order.
  std::string block_str;
};

using VarPtmFilterPtr = std::shared_ptr<VarPtmFilter>;

} /* namespace toppic */

#endif
