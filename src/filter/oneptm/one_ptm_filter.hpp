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

#ifndef TOPPIC_FILTER_ONEPTM_ONE_PTM_FILTER_HPP_
#define TOPPIC_FILTER_ONEPTM_ONE_PTM_FILTER_HPP_

#include <vector>

#include "filter/massmatch/filter_spec.hpp"
#include "filter/massmatch/mass_match.hpp"
#include "filter/massmatch/mass_match_util.hpp"
#include "filter/mng/one_ptm_filter_mng.hpp"
#include "ms/spec/prm_ms.hpp"
#include "prsm/simple_prsm.hpp"
#include "seq/proteoform.hpp"

namespace toppic {

class OnePtmFilter {
 public:
  OnePtmFilter(const ProteoformPtrVec& proteo_ptrs,
               const OnePtmFilterMngPtr& mng_ptr, std::string block_str);
  void computeBestMatch(const PrmMsPtrVec& prm_ms_ptr_vec,
                        const PrmMsPtrVec& srm_ms_ptr_vec);

  // the block independent part of a spectrum set, computed once per stage
  static FilterSpecPtr geneFilterSpec(const PrmMsPtrVec& prm_ms_ptr_vec,
                                     const PrmMsPtrVec& srm_ms_ptr_vec,
                                     const OnePtmFilterMngPtr& mng_ptr);

  void computeBestMatch(const FilterSpecPtr& spec_ptr);

  SimplePrsmPtrVec getCompMatchPtrs() { return comp_match_ptrs_; }
  SimplePrsmPtrVec getPrefMatchPtrs() { return pref_match_ptrs_; }
  SimplePrsmPtrVec getSuffMatchPtrs() { return suff_match_ptrs_; }
  SimplePrsmPtrVec getInternalMatchPtrs() { return internal_match_ptrs_; }

 private:
  OnePtmFilterMngPtr mng_ptr_;
  ProteoformPtrVec proteo_ptrs_;

  MassMatchPtr diag_index_ptr_;
  MassMatchPtr rev_diag_index_ptr_;
  MassMatchPtr term_index_ptr_;
  MassMatchPtr rev_term_index_ptr_;

  // mass_match_util::compNonincreasingShifts of the two N-terminal side
  // indexes, computed once for the block
  std::vector<char> term_nonincr_shifts_;
  std::vector<char> diag_nonincr_shifts_;
  // the sequence ranges of the index rows
  mass_match_util::IndexRowRanges row_ranges_;

  // Score arrays of the four indexes and their per-proteoform maxima. They
  // are reused by the spectra of the database block, not reallocated for each
  // spectrum.
  std::vector<short> term_scores_;
  std::vector<short> diag_scores_;
  std::vector<short> rev_term_scores_;
  std::vector<short> rev_diag_scores_;
  std::vector<int> term_max_;
  std::vector<int> diag_max_;
  std::vector<int> rev_term_max_;
  std::vector<int> rev_diag_max_;

  SimplePrsmPtrVec comp_match_ptrs_;
  SimplePrsmPtrVec pref_match_ptrs_;
  SimplePrsmPtrVec suff_match_ptrs_;
  SimplePrsmPtrVec internal_match_ptrs_;
};

using OnePtmFilterPtr = std::shared_ptr<OnePtmFilter>;
} /* namespace toppic */

#endif
