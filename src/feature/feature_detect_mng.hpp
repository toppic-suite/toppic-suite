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


#ifndef PROT_FEATURE_DETECT_MNG_HPP_
#define PROT_FEATURE_DETECT_MNG_HPP_

#include <vector>
#include "spec/peak_tolerance.hpp"

namespace prot {

class FeatureDetectMng {
 public:
  FeatureDetectMng();

  std::vector<double> getExtMasses(double mass);

  PeakTolerancePtr peak_tolerance_ptr_;

  std::vector<double> ext_offsets_;

  double extend_min_mass_ = 5000;

  int intv_width_ = 500;

};

typedef std::shared_ptr<FeatureDetectMng> FeatureDetectMngPtr;

} /* namespace */

#endif 
