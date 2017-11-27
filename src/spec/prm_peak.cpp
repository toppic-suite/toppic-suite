//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#include <algorithm>
#include <iostream>

#include "base/logger.hpp"
#include "base/support_peak_type_base.hpp"
#include "base/base_data.hpp"
#include "spec/prm_peak.hpp"

namespace prot {

void PrmPeak::addNghbEdge(DeconvPeakPtr deconv_peak_ptr, double offset,
                          SPTypePtr peak_type_ptr, double score) {
  score_ += score;
  SupportPeakPtr support_peak_ptr
      = std::make_shared<SupportPeak>(deconv_peak_ptr, offset, score, peak_type_ptr);
  neighbor_list_.push_back(support_peak_ptr);
}

RmBreakTypePtr PrmPeak::getBreakType() {
  RmBreakTypePtr bt_ptr = RmBreakType::NONE;
  for (size_t i = 0; i < neighbor_list_.size(); i++) {
    if (neighbor_list_[i]->getPeakTypePtr() == SPTypeBase::getSPTypePtr_N_TERM()) {
      if (bt_ptr == RmBreakType::NONE) {
        bt_ptr = RmBreakType::N_TERM;
      } else if (bt_ptr == RmBreakType::C_TERM) {
        bt_ptr = RmBreakType::BOTH;
      }
    } else {
      if (bt_ptr == RmBreakType::NONE) {
        bt_ptr = RmBreakType::C_TERM;
      } else if (bt_ptr == RmBreakType::N_TERM) {
        bt_ptr = RmBreakType::BOTH;
      }
    }
  }
  return bt_ptr;
}

} /* namespace prot */
