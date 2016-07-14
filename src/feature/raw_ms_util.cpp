#ifndef PROT_SPEC_PEAK_LIST_HPP_
#define PROT_SPEC_PEAK_LIST_HPP_

#include <cstddef>
#include <vector>
#include <cmath>

#include "base/logger.hpp"
#include "feature/raw_ms_util.hpp"

namespace prot {

void RawMsUtil::sortOnPos(PeakPtrVec &ptr_list) {
  for (size_t i = 0; i < ptr_list.size(); i++) {
    for (size_t j = i + 1; j < ptr_list.size(); j++) {
      if (ptr_list[i]->getPosition() > ptr_list[j]->getPosition()) {
        PeakPtr tmp = ptr_list[i];
        ptr_list[i] = ptr_list[j];
        ptr_list[j] = tmp;
      }
    }
  }
}

double RawMsUtil::findMaxPos(PeakPtrVec &ptr_list) {
  sortOnPos(ptr_list);
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
