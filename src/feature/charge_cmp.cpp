// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include "feature/charge_cmp.hpp" 

namespace prot {



// get common peaks in two peak lists 
std::vector<int> getCommonPeak(RealEnvPtr env_a, RealEnvPtr env_b) {
  std::vector<int> list_a = env_a->getPeakIdxList();
  std::vector<int> list_b = env_b->getPeakIdxList();
  int len = list_a.size();
  std::vector<int> common_peaks(len, -1);
  for (int i = 0; i < len; i++) {
    for (size_t j = 0; j < list_b.size(); j++) {
      if (env_a->isExist(i) && list_a[i] == list_b[j]) {
        common_peaks[i] = list_a[i];
      }
    }
  }
  return common_peaks;
}

// count the number of common peaks. 
int cntCommonPeak(std::vector<int> &list_a) {
  int cnt = 0;
  for (size_t i = 0; i < list_a.size(); i++) {
    if (list_a[i] >= 0) {
      cnt++;
    }
  }
  return cnt;
}

// check if the second charge state is better 
bool isSecondBetter(PeakPtrVec &peak_list, MatchEnvPtr a,
                    MatchEnvPtr  b, double tolerance) {
  // get common peak 
  std::vector<int> common_peaks = getCommonPeak(a->getRealEnvPtr(), b->getRealEnvPtr());
  int common_num = cntCommonPeak(common_peaks);
  if (common_num <= 6) {
    return false;
  }
  /* get distance list */
  std::vector<double> dist;
  for (size_t i = 0; i < common_peaks.size(); i++) {
    for (size_t j = i + 1; j < common_peaks.size(); j++) {
      if (common_peaks[i] >= 0 && common_peaks[j] >= 0) {
        dist.push_back(
            (peak_list[common_peaks[j]]->getPosition() - peak_list[common_peaks[i]]->getPosition()) / (j - i));
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

int ChargeCmp::comp(PeakPtrVec &peak_list, MatchEnvPtr a, MatchEnvPtr b, double tolerance) {
  if (isSecondBetter(peak_list, a, b, tolerance)) {
    return -1;
  }
  if (isSecondBetter(peak_list, b, a, tolerance)) {
    return 1;
  }
  return 0;
}

}
