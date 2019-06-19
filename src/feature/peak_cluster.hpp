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

#include "spec/raw_ms.hpp"
#include "deconv/env/envelope.hpp"
#include "deconv/env/real_env.hpp"
#include "deconv/env/match_env.hpp"

#include "feature/savitzky_golay.hpp"
#include "feature/frac_feature.hpp"

namespace toppic {

class PeakCluster:public FracFeature {
 public:
  PeakCluster(MatchEnvPtr match_env);

  void addEnvelopes(int min_charge, int max_charge, 
                    int min_ms1_id, int max_ms1_id, 
                    int min_scan_num, int max_scan_num,
                    RealEnvPtrVec envs);

  void clearScores();

  void updateScore(RawMsPtrVec spec_list, bool check_pvalue);

 private:
  // promex feature variables
  int min_ms1_id_;
  int max_ms1_id_;
  double rep_mass_;
  int rep_charge_;
  int rep_ms1_id_;
  double rep_mz_;
  double score_ = 0.0;

  double apex_time_;
  int apex_scan_;
  double apex_intensity_;
  double boundary_intensity_;

  RealEnvPtrVec2D real_envs_;  
  EnvelopePtr theo_env_;

  std::vector<double> rep_summed_intensities_;

  // one for even charge, one for odd charge
  std::vector<int> best_charges_;
  std::vector<double> env_dist_scores_;
  std::vector<double> env_corr_scores_;
  std::vector<double> env_inte_scores_;

  std::vector<double> inte_distr_;
  std::vector<double> best_corr_scores_;
  std::vector<double> best_dist_scores_;
  std::vector<double> best_inte_scores_;
  std::vector<double> xic_corr_between_best_charges_;

  static int even_charge_idx_;
  static int odd_charge_idx_;
  static double win_size_;

  bool init_score_;

  SavitzkyGolayPtr smoother_;

  // do not know the meaning
  int flag_;

};

}
#endif
