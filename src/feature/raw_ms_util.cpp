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


#ifndef PROT_SPEC_PEAK_LIST_HPP_
#define PROT_SPEC_PEAK_LIST_HPP_

#include <cstddef>
#include <vector>
#include <cmath>

#include "base/logger.hpp"
#include "feature/raw_ms_util.hpp"

namespace toppic {

namespace raw_ms_util {

double findMaxPos(const PeakPtrVec &ptr_list) {
  return ptr_list[ptr_list.size() -1]->getPosition();
}

int searchPos(const PeakPtrVec &ptr_list, double pos) {
  if (ptr_list.size() == 0) {
    return -1;
  }
  if (ptr_list.size() == 1) {
    return 0;
  }

  int min = 0;
  int max = ptr_list.size() - 1;
  while (max > min + 1) {
    int mid = (max+min)/2;
    if(ptr_list[mid]->getPosition() ==  pos){
      max = mid;
      min = mid;
    } else if(ptr_list[mid]->getPosition() < pos ){
      min = mid;
    } else if(ptr_list[mid]->getPosition() > pos ){
      max = mid;
    }
  }

  double a = std::abs(ptr_list[max]->getPosition() - pos);
  double b = std::abs(ptr_list[min]->getPosition() - pos);
  if (a < b) {
    return max;
  } else {
    return min;
  }
}

// Finds the nearest peak for a specific position with error tolerance.
int getNearPeakIdx(const PeakPtrVec &ptr_list, double pos, double tolerance) {
  // find the peak nearest to pos 
  int idx = searchPos(ptr_list, pos);
  //LOG_DEBUG("pos " << pos << " peak pos " << ptr_list[idx]->getPosition());
  if (idx < 0 || std::abs(ptr_list[idx]->getPosition() - pos) > tolerance) {
    return -1;
  } else {
    return idx;
  }
}

// Removes a list of peaks.
PeakPtrVec rmPeaks(const PeakPtrVec &ptr_list, std::vector<bool> &keep) {
  PeakPtrVec new_list;
  for (size_t i = 0; i < ptr_list.size(); i++) {
    if (keep[i]) {
      new_list.push_back(ptr_list[i]);
    }
  }
  return new_list;
}

}  // namespace raw_ms_util

}  // namespace toppic

#endif
