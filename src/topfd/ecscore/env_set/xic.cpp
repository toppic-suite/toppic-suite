//Copyright (c) 2014 - 2025, The Trustees of Indiana University.
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

#include <numeric>
#include "topfd/ecscore/env_set/xic.hpp"

namespace toppic {

Xic::Xic(std::vector<double> inte_ratio_list,
         std::vector<double> top_three_inte_list,
         std::vector<double> all_peak_inte_list) {
  inte_ratio_list_ = inte_ratio_list;
  top_three_inte_list_ = top_three_inte_list;
  top_three_inte_sum_ = std::accumulate(top_three_inte_list_.begin(), 
                                        top_three_inte_list_.end(), 
                                        0.0); 
  all_peak_inte_list_ = all_peak_inte_list;
  all_peak_inte_sum_ = std::accumulate(all_peak_inte_list_.begin(),
                                       all_peak_inte_list_.end(),
                                       0.0);
  moving_avg(moving_avg_win_);
}

void Xic::moving_avg(int size) {
  std::vector<double> data{0};
  data.insert(data.end(), top_three_inte_list_.begin(), top_three_inte_list_.end());
  int num_spec = data.size();
  double sum = 0.0;
  int cnt = 0;
  for (int i = 0; i < num_spec; i++) {
    sum += data[i];
    cnt++;
    if (cnt >= size) {
      smoothed_inte_list_.push_back((sum / (double) size));
      sum -= data[cnt - size];
    }
  }
}

}
