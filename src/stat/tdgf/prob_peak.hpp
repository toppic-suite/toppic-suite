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

#ifndef TOPPIC_STAT_TDGF_PROB_PEAK_HPP_
#define TOPPIC_STAT_TDGF_PROB_PEAK_HPP_

#include "ms/spec/prm_peak.hpp"
#include "ms/spec/base_peak_type.hpp"

namespace toppic {

class ProbPeak {
 public:
  ProbPeak(PrmPeakPtr peak_ptr, int spectrum_id, int height, 
           bool strict, double convert_ratio);
  int mass_;
  int tolerance_;
  BasePeakTypePtr base_type_ptr_;
  int spectrum_id_;
  int mass_bgn_;
  int mass_end_;
  int table_bgn_;
  int table_end_;
};

}

#endif
