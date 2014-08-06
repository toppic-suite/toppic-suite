#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

namespace prot {

template <class T>
T findHighestPeak(std::vector<T> ptr_list) {
  int idx = findHighestPeakIdx(ptr_list);
  if (idx >= 0) {
    return ptr_list[idx];
  } else {
    return T(nullptr);
  }
}

template <class T>
int findHighestPeakIdx(const std::vector<T> &ptr_list) {
  int idx = -1;
  for (size_t i = 0; i < ptr_list.size(); i++) {
    if (idx < 0 
        || ptr_list[i]->getIntensity() > ptr_list[idx]->getIntensity()) {
      idx = i;
    }
  }
  return idx;
}

template <class T>
void sortOnPos(std::vector<T> &ptr_list) {
  for (size_t i = 0; i < ptr_list.size(); i++) {
    for (size_t j = i + 1; j < ptr_list.size(); j++) {
      if (ptr_list[i]->getPosition() > ptr_list[j]->getPosition()) {
        T tmp = ptr_list[i];
        ptr_list[i] = ptr_list[j];
        ptr_list[j] = tmp;
      }
    }
  }
}

template <class T>
double findMaxPos(const std::vector<T> &ptr_list) {
  sortOnPos(ptr_list);
  return ptr_list[ptr_list.size() -1]->getPosition();
}

template <class T>
void sortOnIntensity(std::vector<T> &ptr_list) {
  for (size_t i = 0; i < ptr_list.size(); i++) {
    for (size_t j = i + 1; j < ptr_list.size(); j++) {
      if (ptr_list[i]->getIntensity() < ptr_list[j]->getIntensity()) {
        T tmp = ptr_list[i];
        ptr_list[i] = ptr_list[j];
        ptr_list[j] = tmp;
      }
    }
  }
}

template <class T>
int searchPos(const std::vector<T> &ptr_list, double pos) {
  int idx = -1;
  int min = 0;
  int max = ptr_list.size() - 1;
  while (max >= min && pos == -1) {
    int mid = (max+min)/2;
    if(ptr_list[mid]->getPosition() ==  pos){
      idx = mid;
    }else if(ptr_list[mid]->getPostion() < pos){
      min = mid +1;
    }else if(ptr_list[mid] > pos){
      max = mid -1;
    }
  }
  if (idx >= 0) {
    return idx;
  }
  else {
    if (ptr_list.size() == 0) {
      return -1;
    }
    else if (ptr_list.size() == 1) {
      return 0;
    }
    else {
      double a = std::abs(ptr_list[max] - pos);
      double b = std::abs(ptr_list[min] - pos);
      if (a < b) {
        return max;
      }
      else {
        return min;
      }
    }
  }
}


/**
 * Finds the highest peak for a specific position with error tolerance.
 */
template <class T>
int getHighPeakIdx(const std::vector<T> &ptr_list, double pos, double tolerance) {
  int idx = searchPos(ptr_list, pos);
  int best_idx = -1;
  /* extend to left */
  int i = idx - 1;
  while (i >= 0 
         && std::abs(ptr_list[i].getPosition() - pos) <= tolerance) {
    if (best_idx < 0
        || ptr_list[i].getIntensity() > ptr_list[best_idx].getIntensity()) {
      best_idx = i;
    }
    i--;
  }
  /* extend to right */
  i = idx;
  while (i < ptr_list.size() 
         && std::abs(ptr_list[i].getPosition() - pos) <= tolerance) {
    if (best_idx < 0
        || ptr_list[i].getIntensity() > ptr_list[best_idx].getIntensity()) {
      best_idx = i;
    }
    i++;
  }
  return best_idx;
}

/**
 * Finds the nearest peak for a specific position with error tolerance.
 */
template <class T>
int getNearPeakIdx(const std::vector<T> &ptr_list, double pos, double tolerance) {
  /* find the peak nearest to pos */
  int idx = searchPos(ptr_list, pos);
  if (idx < 0 || std::abs(ptr_list[idx]->getPosition - pos) > tolerance) {
    return -1;
  }
  else {
    return idx;
  }
}

/**
 * Removes all peaks in [center-interval, center+interval]
 */
template <class T>
std::vector<T> rmPeaks(std::vector<T> &ptr_list, double center, 
                       double interval) {
  std::vector<T> new_list;
  for (size_t i = 0; i < ptr_list.size(); i++) {
    if (std::abs(ptr_list[i].getPosition() - center) > interval) {
      new_list.push_back(ptr_list[i]);
    }
  }
  return new_list;
}

template <class T>
std::vector<T> rmZeroPeaks(std::vector<T> &ptr_list) {
  std::vector<T> new_list;
  for (size_t i = 0; i < ptr_list.size(); i++) {
    if (ptr_list[i].getIntensity() != 0) {
      new_list.push_back(ptr_list[i]);
    }
  }
  return new_list;
}

/**
 * Removes a list of peaks.
 */
template <class T>
std::vector<T> rmPeaks(std::vector<T> &ptr_list, std::vector<bool> keep) {
  std::vector<T> new_list;
  for (size_t i = 0; i < ptr_list.size(); i++) {
    if (keep[i]) {
      new_list.push_back(ptr_list[i]);
    }
  }
  return new_list;
}

/**
 * Removes precursor mass. In top down data, one peak may be split to two
 * peaks. try to combine them together.
 */
template <class T>
std::vector<T> rmClosePeaks(std::vector<T> &ptr_list, double tolerance) {
  sortOnPos(ptr_list);
  for (size_t i = 0; i < ptr_list.size() - 1; i++) {
    T pA = ptr_list[i];
    T pB = ptr_list[i + 1];
    if (std::abs(pA->getPosition() - pB->getPosition()) <= tolerance) {
      double a = pA->getIntensity();
      double b = pB->getIntensity();
      if (a > b) {
        pA->setIntensity(a + b);
        pB->setIntensity(0);
      } else {
        pA->setIntensity(0);
        pB->setIntensity(a + b);
      }
    }
  }
  return rmZeroPeaks(ptr_list);
}

}

