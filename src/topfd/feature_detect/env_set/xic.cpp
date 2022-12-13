//Copyright (c) 2014 - 2022, The Trustees of Indiana University.
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


#include "xic.hpp"

namespace toppic {
  Xic::Xic() {
    start_spec_id_ = -1;
    base_spec_id_ = -1;
  }
  Xic::Xic(int start_spec_id, int base_spec_id, std::vector<double> &inte_list,
                   std::vector<double> &env_inte_list) {
    start_spec_id_ = start_spec_id;
    base_spec_id_ = base_spec_id;
    inte_list_ = inte_list;
    env_inte_list_ = env_inte_list;
    moving_avg(2);
  }

  Xic::Xic(int start_spec_id, int base_spec_id, std::vector<double> &inte_list,
                   std::vector<double> &smoothed_inte_list, std::vector<double> &env_inte_list) {
    start_spec_id_ = start_spec_id;
    base_spec_id_ = base_spec_id;
    inte_list_ = inte_list;
    smoothed_inte_list_ = smoothed_inte_list;
    env_inte_list_ = env_inte_list;
  }

  Xic::Xic(const Xic &x) {
    start_spec_id_ = x.start_spec_id_;
    base_spec_id_ = x.base_spec_id_;
    for (auto inte: x.inte_list_)
      inte_list_.push_back(inte);
    for (auto smoothed_inte: x.smoothed_inte_list_)
      smoothed_inte_list_.push_back(smoothed_inte);
    for (auto env_inte: x.env_inte_list_)
      env_inte_list_.push_back(env_inte);
  }

  void Xic::moving_avg(int size) {
    std::vector<double> data = inte_list_;
    std::vector<double> left_padding(1, 0);
    data.insert(data.begin(), left_padding.begin(), left_padding.end());
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

  bool Xic::isEmpty() {
    if (start_spec_id_ == -1 and base_spec_id_ == -1 and inte_list_.empty() and smoothed_inte_list_.empty())
      return true;
    return false;
  }
}