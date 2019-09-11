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


#ifndef TOPPIC_FEATURE_PEAK_CLUSTER_HPP_
#define TOPPIC_FEATURE_PEAK_CLUSTER_HPP_

#include <memory>
#include <vector>

#include "spec/raw_ms.hpp"
#include "deconv/env/envelope.hpp"
#include "deconv/env/real_env.hpp"

#include "feature/frac_feature.hpp"
#include "feature/savitzky_golay.hpp"

namespace toppic {

class PeakCluster:public FracFeature {
 public:
  PeakCluster(EnvelopePtr theo_env);

  void addEnvelopes(FracFeaturePtr feature_ptr, 
                    RealEnvPtrVec envs);

  void clearScores();

  void updateScore(PeakPtrVec2D &raw_peaks, bool check_pvalue);

  double getInteDistr(int i) {return inte_distr_[i];}

  double getSumDistScore(int i) {return sum_dist_scores_[i];}
  double getSumCorrScore(int i) {return sum_corr_scores_[i];}
  double getSumInteScore(int i) {return sum_inte_scores_[i];}

  double getBestDistScore(int i) {return best_dist_scores_[i];}
  double getBestCorrScore(int i) {return best_corr_scores_[i];}
  double getBestInteScore(int i) {return best_inte_scores_[i];}

  double getXicCorrBetweenCharges(int i) {return xic_corr_between_best_charges_[i];}

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
  std::vector<double> sum_dist_scores_;
  std::vector<double> sum_corr_scores_;
  std::vector<double> sum_inte_scores_;

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

typedef std::shared_ptr<PeakCluster> PeakClusterPtr;

}
#endif
