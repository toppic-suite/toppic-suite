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

#include "exp_envelope.hpp"

namespace toppic {
  ExpEnvelope::ExpEnvelope() {
    spec_id_ = -1;
    peak_list_.empty();
  }

  ExpEnvelope::ExpEnvelope(int spec_id, std::vector<ExpPeak> peak_list) {
    spec_id_ = spec_id;
    peak_list_ = peak_list;
  }

  int ExpEnvelope::get_match_peak_num(int base_idx) {
    int total_peaks = peak_list_.size();
    int start_idx = std::max(base_idx - 1, 0);
    int end_idx = std::min(base_idx + 1, total_peaks - 1);
    int num = 0;
    for (int i = start_idx; i < end_idx + 1; i++) {
      ExpPeak p = peak_list_[i];
      if (!p.isEmpty())
        num = num + 1;
    }
    return num;
  }

  std::vector<double> ExpEnvelope::get_inte_list() {
    std::vector<double> inte_list;
    for (auto p: peak_list_) {
      if (!p.isEmpty())
        inte_list.push_back(p.getInte());
      else
        inte_list.push_back(0);
    }
    return inte_list;
  }

  std::vector<double> ExpEnvelope::get_pos_list() {
    std::vector<double> pos_list;
    for (auto p: peak_list_) {
      if (!p.isEmpty())
        pos_list.push_back(p.getPos());
      else
        pos_list.push_back(0);
    }
    return pos_list;
  }

  std::vector<double> ExpEnvelope::get_non_empty_pos_list() {
    std::vector<double> pos_list;
    for (auto p: peak_list_) {
      if (!p.isEmpty())
        pos_list.push_back(p.getPos());
    }
    return pos_list;
  }

  void ExpEnvelope::setExpEnvList(const std::vector<ExpPeak> &peak_list) {
    for (size_t i = 0; i < peak_list.size(); i++)
      peak_list_[i] = peak_list[i];
  }

  bool ExpEnvelope::isEmpty() {
    if (spec_id_ == -1 && peak_list_.empty())
      return true;
    return false;
  }
}
