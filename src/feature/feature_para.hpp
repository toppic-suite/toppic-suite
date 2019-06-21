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


#ifndef TOPPIC_FEATURE_FEATURE_PARA_HPP_
#define TOPPIC_FEATURE_FEATURE_PARA_HPP_

#include <vector>
#include "spec/peak_tolerance.hpp"
#include "feature/peak_cluster_score.hpp"

namespace toppic {

class FeaturePara {
 public:
  FeaturePara(int frac_id, const std::string &file_name, 
              std::string &resourse_dir); 

  std::vector<double> getExtendMasses(double mass);

  std::vector<double> getExtendOffsets() {return extend_offsets_;}

  std::vector<double> getSearchMasses(double mass);

  std::vector<double> getSearchOffsets() {return search_offsets_;}

  PeakTolerancePtr peak_tolerance_ptr_;

  // extend_offset is -2, -1, 0, 1, 2. We use it to search matched masses 
  // when a reference mass is given to find a feature. 
  std::vector<double> extend_offsets_;

  // search offset is  -3, -2, -1, 0, 1, 2, 3. We use it to search the precursor 
  // mass of an MS/MS spectrum to find a matched feature.
  std::vector<double> search_offsets_;

  PeakClusterScorePtr peak_cluster_score_ptr_;

  double extend_min_mass_ = 5000;

  int intv_width_ = 500;

  int feature_num_ = 20000;

  int frac_id_;

  std::string file_name_;

};

typedef std::shared_ptr<FeaturePara> FeatureParaPtr;

} /* namespace */

#endif 
