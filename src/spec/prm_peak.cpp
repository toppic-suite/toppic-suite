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


#include <algorithm>
#include <iostream>

#include "base/logger.hpp"
#include "base/support_peak_type_base.hpp"
#include "base/base_data.hpp"
#include "spec/prm_peak.hpp"

namespace prot {

void PrmPeak::addNghbEdge(DeconvPeakPtr deconv_peak_ptr,double offset,
                          SPTypePtr peak_type_ptr,double score){
  score_ += score;
  SupportPeakPtr support_peak_ptr(
      new SupportPeak(deconv_peak_ptr,offset,score, peak_type_ptr));
  neighbor_list_.push_back(support_peak_ptr);
}

RmBreakTypePtr PrmPeak::getBreakType() {
  RmBreakTypePtr bt_ptr = RmBreakType::NONE;
  for(size_t i=0; i<neighbor_list_.size(); i++){
    if(neighbor_list_[i]->getPeakTypePtr() == 
       SPTypeBase::getSPTypePtr_N_TERM()) {
      if(bt_ptr == RmBreakType::NONE) {
        bt_ptr = RmBreakType::N_TERM;
      }
      else if(bt_ptr == RmBreakType::C_TERM){
        bt_ptr = RmBreakType::BOTH;
      }
    }
    else{
      if(bt_ptr == RmBreakType::NONE){
        bt_ptr = RmBreakType::C_TERM;
      }
      else if(bt_ptr == RmBreakType::N_TERM){
        bt_ptr = RmBreakType::BOTH;
      }
    }
  }
  return bt_ptr;
}

} /* namespace prot */
