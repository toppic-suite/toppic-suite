#ifndef PROT_FEATURE_RAW_MS_UTIL_HPP_
#define PROT_FEATURE_RAW_MS_UTIL_HPP_

#include <cstddef>
#include <vector>
#include <cmath>

#include "spec/peak.hpp"

namespace prot {

class RawMsUtil {
 public: 
  static void sortOnPos(PeakPtrVec &ptr_list);

  static double findMaxPos(PeakPtrVec &ptr_list);

  static int searchPos(PeakPtrVec &ptr_list, double pos);

  static int getNearPeakIdx(PeakPtrVec  &ptr_list, double pos, double tolerance);

  static PeakPtrVec rmPeaks(PeakPtrVec &ptr_list, std::vector<bool> &keep);

};

}

#endif
