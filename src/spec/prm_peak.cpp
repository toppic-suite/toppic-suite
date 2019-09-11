//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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


#include "common/util/logger.hpp"
#include "common/base/support_peak_type_base.hpp"
#include "spec/prm_peak.hpp"

namespace toppic {

PrmPeak::PrmPeak(int spec_id, DeconvPeakPtr base_peak_ptr,
          BasePeakTypePtr base_type,
          double mono_mass, double score,
          double strict_tolerance,
          double n_strict_c_relax_tolerance,
          double n_relax_c_strict_tolerance):
      Peak(mono_mass, base_peak_ptr->getIntensity()),
      spec_id_(spec_id),
      base_peak_ptr_(base_peak_ptr),
      base_type_(base_type),
      mono_mass_(mono_mass),
      score_(score),
      strict_tolerance_(strict_tolerance),
      n_strict_c_relax_tolerance_(n_strict_c_relax_tolerance),
      n_relax_c_strict_tolerance_(n_relax_c_strict_tolerance) {}

void PrmPeak::setMonoMass(double m) {
  mono_mass_ = m;
  setPosition(m);
}

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

} /* namespace toppic */
