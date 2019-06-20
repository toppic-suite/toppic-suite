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


#ifndef TOPPIC_FEATURE_PEAK_CLUSTER_SCORE_HPP_
#define TOPPIC_FEATURE_PEAK_CLUSTER_SCORE_HPP_

#include <memory>
#include <vector>

#include "feature/peak_cluster.hpp"

namespace toppic {

class PeakClusterScore {
 public:
  PeakClusterScore(std::string &dir, double thresh); 

  double getScore(PeakClusterPtr pc); 

 private:
  std::vector<std::vector<double>> dist_score_table_;
  std::vector<std::vector<double>> corr_score_table_;
  std::vector<std::vector<double>> inte_score_table_;

  std::vector<std::vector<double>> dist_sum_score_table_;
  std::vector<std::vector<double>> corr_sum_score_table_;
  std::vector<std::vector<double>> inte_sum_score_table_;

  std::vector<std::vector<double>> xic_score_table_1_;
  std::vector<std::vector<double>> xic_score_table_2_;

  std::vector<std::vector<double>> inte_distr_score_table_;

  std::vector<double> mass_bins_;

  // row_num, column number are fixed
  int row_num_ = 30;
  int col_num_  = 1001;

  double score_thresh_;
};

}
#endif
