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

#include <cmath>

#include "stat/tdgf/prob_peak.hpp"

namespace toppic {

ProbPeak::ProbPeak(PrmPeakPtr peak_ptr, int spectrum_id, int height, 
                   bool strict, double convert_ratio) {

  mass_ = static_cast<int>(std::round(peak_ptr->getMonoMass() * convert_ratio));
  // we use NStrict and CRelax tolerance
  if (strict) {
    tolerance_ = std::ceil(peak_ptr->getStrictTolerance() * convert_ratio);
  } else {
    tolerance_ = std::ceil(peak_ptr->getNStrictCRelaxTolerance() * convert_ratio);
  }
  base_type_ptr_ = peak_ptr->getBaseTypePtr();
  spectrum_id_ = spectrum_id;
  mass_bgn_ = mass_ - tolerance_;
  if (mass_bgn_ < 0) {
    mass_bgn_ = 0;
  }
  mass_end_ = mass_ + tolerance_;
  table_bgn_ = mass_bgn_ * height;
  table_end_ = mass_end_ * height + (height - 1);
}

}
