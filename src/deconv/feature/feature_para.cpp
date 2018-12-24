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


#include "base/mass_constant.hpp"
#include "deconv/feature/feature_para.hpp"

namespace toppic {

FeaturePara::FeaturePara() {
  double ppo = 0.000015;
  bool use_min_tolerance = true;
  double min_tolerance = 0.01;
  peak_tolerance_ptr_
      = std::make_shared<PeakTolerance>(ppo, use_min_tolerance, min_tolerance);

  // extend sp parameter 
  double IM = mass_constant::getIsotopeMass();
  // the set of offsets used to expand the monoisotopic mass list 
  std::vector<double> offsets {{0, -IM, IM, -2 * IM, 2 * IM}};
  ext_offsets_ = offsets;
}

std::vector<double> FeaturePara::getExtMasses(double mass) {
  std::vector<double> result;
  if (mass < extend_min_mass_) {
    result.push_back(mass);
  }
  else {
    for (size_t i = 0; i < ext_offsets_.size(); i++) {
      double new_mass = mass + ext_offsets_[i];
      result.push_back(new_mass);
    }
  }
  return result;
}

} /* namespace */

