// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#ifndef PROT_PTM_SEARCH_SLOW_FILTER_HPP_
#define PROT_PTM_SEARCH_SLOW_FILTER_HPP_

#include <memory>

#include "spec/spectrum_set.hpp"
#include "prsm/simple_prsm.hpp"
#include "oneptmsearch/ptm_search_mng.hpp"
#include "ptmsearch/ptm_slow_match.hpp"
#include "ptmsearch/comp_shift_low_mem.hpp"

namespace prot {

class PtmSearchSlowFilter {
 public:
  PtmSearchSlowFilter(SpectrumSetPtr spectrum_set_ptr,
                      SimplePrsmPtrVec simple_prsm_ptrs,
                      CompShiftLowMem comp_shift,
                      PtmSearchMngPtr mng_ptr);
  PrsmPtrVec getPrsms(int shift_num, AlignTypePtr type_ptr);

 private:
  PtmSlowMatchPtrVec complete_prefix_slow_match_ptrs_;
  PtmSlowMatchPtrVec suffix_internal_slow_match_ptrs_;
  PrsmPtrVec2D complete_prsm_2d_ptrs_;
  PrsmPtrVec2D prefix_prsm_2d_ptrs_;
  PrsmPtrVec2D suffix_prsm_2d_ptrs_;
  PrsmPtrVec2D internal_prsm_2d_ptrs_;
};

typedef std::shared_ptr<PtmSearchSlowFilter> PtmSearchSlowFilterPtr;

} /* namespace prot */

#endif /* PTM_SEARCH_SLOW_FILTER_HPP_ */
