//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#include <cmath>

#include "ms/env/charge_cmp.hpp" 

namespace toppic {

namespace charge_cmp {

// get shared peaks in two peak lists 
std::vector<int> getSharedPeaks(ExpEnvPtr env_a, ExpEnvPtr env_b) {
  int len = env_a->getPeakNum();
  std::vector<int> shared_peaks(len, -1);
  for (int i = 0; i < len; i++) {
    int a_idx = env_a->getPeakIdx(i);
    for (int j = 0; j < env_b->getPeakNum(); j++) {
      int b_idx = env_b->getPeakIdx(j);
      if (env_a->isExist(i) && a_idx == b_idx) {
        shared_peaks[i] = a_idx;
      }
    }
  }
  return shared_peaks;
}

// count the number of shared peaks. 
int cntSharedPeaks(std::vector<int> &list_a) {
  int cnt = 0;
  for (size_t i = 0; i < list_a.size(); i++) {
    if (list_a[i] >= 0) {
      cnt++;
    }
  }
  return cnt;
}

// check if the second charge state is better 
bool isSecondBetter(const PeakPtrVec &peak_list, MatchEnvPtr a,
                    MatchEnvPtr  b, double tolerance) {
  // get shared peak 
  std::vector<int> s_peaks = getSharedPeaks(a->getRealEnvPtr(), 
                                                 b->getRealEnvPtr());
  int s_num = cntSharedPeaks(s_peaks);
  if (s_num <= 6) {
    return false;
  }
  /* get distance list */
  std::vector<double> dist;
  double d; 
  for (size_t i = 0; i < s_peaks.size(); i++) {
    for (size_t j = i + 1; j < s_peaks.size(); j++) {
      if (s_peaks[i] >= 0 && s_peaks[j] >= 0) {
        d = (peak_list[s_peaks[j]]->getPosition() - peak_list[s_peaks[i]]->getPosition()) / (j - i);
        dist.push_back(d);
      }
    }
  }
  // check if the distance is near to 1/(chrg_a) or 1/(chrg+b) 
  int cnt_a = 0;
  int cnt_b = 0;
  double center_a = 1 / a->getTheoEnvPtr()->getCharge();
  double center_b = 1 / b->getTheoEnvPtr()->getCharge();
  for (size_t i = 0; i < dist.size(); i++) {
    double d = dist[i];
    if (std::abs(d - center_a) < std::abs(d - center_b)
        && std::abs(d - center_a) < tolerance) {
      cnt_a++;
    }
    if (std::abs(d - center_b) < std::abs(d - center_a)
        && std::abs(d - center_b) < tolerance) {
      cnt_b++;
    }
  }
  if (cnt_b > cnt_a) {
    return true;
  } else {
    return false;
  }
}

int comp(const PeakPtrVec &peak_list, MatchEnvPtr a, MatchEnvPtr b, double tolerance) {
  if (isSecondBetter(peak_list, a, b, tolerance)) {
    return -1;
  }
  if (isSecondBetter(peak_list, b, a, tolerance)) {
    return 1;
  }
  return 0;
}

}

}
