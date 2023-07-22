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


#include "common/util/logger.hpp"

#include "topfd/ecscore/env/ms_map_env.hpp"

namespace toppic {

MsMapEnv::MsMapEnv(int spec_id, MsMapPeakPtrVec peak_list) {
  spec_id_ = spec_id;
  peak_list_ = peak_list;
}

int MsMapEnv::getTopThreeMatchNum(int ref_idx) {
  int total_peaks = peak_list_.size();
  int start_idx = std::max(ref_idx - 1, 0);
  int end_idx = std::min(ref_idx + 1, total_peaks - 1);
  int num = 0;
  for (int i = start_idx; i < end_idx + 1; i++) {
    MsMapPeakPtr p = peak_list_[i];
    if (p != nullptr) {
      num = num + 1;
    }
  }
  return num;
}

std::vector<double> MsMapEnv::getInteList() {
  std::vector<double> inte_list;
  for (auto p: peak_list_) {
    if (p != nullptr)
      inte_list.push_back(p->getIntensity());
    else
      inte_list.push_back(0);
  }
  return inte_list;
}

std::vector<double> MsMapEnv::getMzList() {
  std::vector<double> pos_list;
  for (auto p: peak_list_) {
    if (p != nullptr)
      pos_list.push_back(p->getPosition());
    else
      pos_list.push_back(0);
  }
  return pos_list;
}

}
