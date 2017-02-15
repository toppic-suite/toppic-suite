// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#ifndef PROT_SPEC_PEAK_LIST_HPP_
#define PROT_SPEC_PEAK_LIST_HPP_

#include <cstddef>
#include <vector>
#include <cmath>

#include "base/logger.hpp"
#include "feature/raw_ms_util.hpp"

namespace prot {

double RawMsUtil::findMaxPos(PeakPtrVec &ptr_list) {
  return ptr_list[ptr_list.size() -1]->getPosition();
}


int RawMsUtil::searchPos(PeakPtrVec &ptr_list, double pos) {
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
    }else if(ptr_list[mid]->getPosition() < pos ){
      min = mid;
    }else if(ptr_list[mid]->getPosition() > pos ){
      max = mid;
    }
  }

  double a = std::abs(ptr_list[max]->getPosition() - pos);
  double b = std::abs(ptr_list[min]->getPosition() - pos);
  if (a < b) {
    return max;
  }
  else {
    return min;
  }
}

// Finds the nearest peak for a specific position with error tolerance.
int RawMsUtil::getNearPeakIdx(PeakPtrVec &ptr_list, double pos, double tolerance) {
  // find the peak nearest to pos 
  int idx = searchPos(ptr_list, pos);
  //LOG_DEBUG("pos " << pos << " peak pos " << ptr_list[idx]->getPosition());
  if (idx < 0 || std::abs(ptr_list[idx]->getPosition() - pos) > tolerance) {
    return -1;
  }
  else {
    return idx;
  }
}

// Removes a list of peaks.
PeakPtrVec RawMsUtil::rmPeaks(PeakPtrVec &ptr_list, std::vector<bool> &keep) {
  PeakPtrVec new_list;
  for (size_t i = 0; i < ptr_list.size(); i++) {
    if (keep[i]) {
      new_list.push_back(ptr_list[i]);
    }
  }
  return new_list;
}


}

#endif
