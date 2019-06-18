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


#ifndef TOPPIC_FEATURE_PEAK_CLUSTER_HPP_
#define TOPPIC_FEATURE_PEAK_CLUSTER_HPP_

#include <memory>
#include <vector>

#include "deconv/env/envelope.hpp"
#include "deconv/env/real_env.hpp"

#include "feature/frac_feature.hpp"

namespace toppic {

class PeakCluster:public FracFeature {
 public:

 private:
  double apex_time_;
  int apex_scan_;
  double apex_intensity_;
  double boundary_intensity_;

  RealEnvPtrVec2D envelopes_;  
  EnvelopePtr theo_env_;

  // one for even charge, one for odd charge
  std::vector<int> best_charges_;
  std::vector<double> rep_summed_envelopes_;
  std::vector<double> env_dist_scores_;
  std::vector<double> env_corr_scores_;
  std::vector<double> env_inte_scores_;

  std::vector<double> inte_distr_;
  std::vector<double> best_corr_scores_;
  std::vector<double> best_dist_scores_;
  std::vector<double> best_inte_scores_;
  std::vector<double> xic_corr_between_best_charges_;
};

}
#endif
