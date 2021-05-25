//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#ifndef TOPPIC_SEARCH_ZERO_PTM_SEARCH_ZERO_PTM_FAST_SEARCH_HPP_
#define TOPPIC_SEARCH_ZERO_PTM_SEARCH_ZERO_PTM_FAST_SEARCH_HPP_

#include "seq/proteoform.hpp"
#include "ms/spec/extend_ms.hpp"
#include "search/zeroptmsearch/zero_ptm_fast_match.hpp"

namespace toppic {

namespace zero_ptm_fast_search {

ZpFastMatchPtrVec filter(ProteoformTypePtr align_type_ptr,
                         const ExtendMsPtrVec &ms_ptr_ptr,
                         const ProteoformPtrVec &proteo_ptrs,
                         int report_num, double ppo);

}

}  // namespace toppic
#endif

