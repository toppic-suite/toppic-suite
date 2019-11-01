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


#include "common/util/file_util.hpp"
#include "common/base/mass_constant.hpp"
#include "feature/feature_para.hpp"

namespace toppic {

FeaturePara::FeaturePara(int frac_id, const std::string &file_name, 
                         const std::string &resource_dir): 
  frac_id_(frac_id),
  file_name_(file_name) {

  double ppo = 0.000015;
  bool use_min_tolerance = true;
  double min_tolerance = 0.01;
  peak_tolerance_ptr_
      = std::make_shared<PeakTolerance>(ppo, use_min_tolerance, min_tolerance);

  // extend sp parameter 
  double IM = mass_constant::getIsotopeMass();
  // the set of offsets used to expand the monoisotopic mass list 
  std::vector<double> offsets_1 {{0, -IM, IM, -2 * IM, 2 * IM, -3*IM, 3*IM}};
  search_offsets_ = offsets_1;
  std::vector<double> offsets_2 {{0, -IM, IM, -2 * IM, 2 * IM}};
  extend_offsets_ = offsets_2;

  //peak_cluster_score
  double threshold = 0;
  std::string dir = resource_dir + file_util::getFileSeparator() + "promex"; 
  peak_cluster_score_ptr_ = std::make_shared<PeakClusterScore>(dir, threshold);
}

std::vector<double> FeaturePara::getExtendMasses(double mass) {
  std::vector<double> result;
  if (mass < extend_min_mass_) {
    result.push_back(mass);
  }
  else {
    for (size_t i = 0; i < extend_offsets_.size(); i++) {
      double new_mass = mass + extend_offsets_[i];
      result.push_back(new_mass);
    }
  }
  return result;
}

std::vector<double> FeaturePara::getSearchMasses(double mass) {
  std::vector<double> result;
  if (mass < extend_min_mass_) {
    result.push_back(mass);
  }
  else {
    for (size_t i = 0; i < search_offsets_.size(); i++) {
      double new_mass = mass + search_offsets_[i];
      result.push_back(new_mass);
    }
  }
  return result;
}

} /* namespace */

