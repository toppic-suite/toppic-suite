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

#ifndef TOPPIC_SEARCH_ZEROPTMSEARCH_ZERO_PTM_FAST_SEARCH_HPP_
#define TOPPIC_SEARCH_ZEROPTMSEARCH_ZERO_PTM_FAST_SEARCH_HPP_

#include <utility>
#include <vector>

#include "ms/spec/extend_ms.hpp"
#include "para/peak_tolerance.hpp"
#include "search/zeroptmsearch/zero_ptm_fast_match.hpp"
#include "seq/proteoform.hpp"

namespace toppic {

namespace zero_ptm_fast_search {

ZpFastMatchPtrVec filter(const ProteoformTypePtr& align_type_ptr,
                         const ExtendMsPtrVec& ms_ptr_ptr,
                         const ProteoformPtrVec& proteo_ptrs, int report_num,
                         const PeakTolerancePtr& prec_tole_ptr);

// The same, with the sequence range of each proteoform given (issue #348):
// ranges[i] is the first and the last residue of the range, 0-based in the
// residues of proteo_ptrs[i]. Only that range is tested against the
// precursor mass and scored; a range of -1 is searched for as filter does.
ZpFastMatchPtrVec filterRanges(const ProteoformTypePtr& align_type_ptr,
                               const ExtendMsPtrVec& ms_ptr_vec,
                               const ProteoformPtrVec& proteo_ptrs,
                               const std::vector<std::pair<int, int>>& ranges,
                               int report_num,
                               const PeakTolerancePtr& prec_tole_ptr);

}

}  // namespace toppic
#endif
