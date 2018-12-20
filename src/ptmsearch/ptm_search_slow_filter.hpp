//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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


#ifndef PROT_PTM_SEARCH_SLOW_FILTER_HPP_
#define PROT_PTM_SEARCH_SLOW_FILTER_HPP_

#include <memory>

#include "spec/spectrum_set.hpp"
#include "prsm/simple_prsm.hpp"
#include "oneptmsearch/ptm_search_mng.hpp"
#include "ptmsearch/ptm_slow_match.hpp"
#include "ptmsearch/comp_shift_low_mem.hpp"

namespace toppic {

class PtmSearchSlowFilter {
 public:
  PtmSearchSlowFilter(SpectrumSetPtr spectrum_set_ptr,
                      SimplePrsmPtrVec simple_prsm_ptrs,
                      CompShiftLowMem comp_shift,
                      PtmSearchMngPtr mng_ptr);
  PrsmPtrVec getPrsms(int shift_num, ProteoformTypePtr type_ptr);

 private:
  PtmSlowMatchPtrVec complete_prefix_slow_match_ptrs_;
  PtmSlowMatchPtrVec suffix_internal_slow_match_ptrs_;
  PrsmPtrVec2D complete_prsm_2d_ptrs_;
  PrsmPtrVec2D prefix_prsm_2d_ptrs_;
  PrsmPtrVec2D suffix_prsm_2d_ptrs_;
  PrsmPtrVec2D internal_prsm_2d_ptrs_;
};

typedef std::shared_ptr<PtmSearchSlowFilter> PtmSearchSlowFilterPtr;

} /* namespace toppic */

#endif /* PTM_SEARCH_SLOW_FILTER_HPP_ */
