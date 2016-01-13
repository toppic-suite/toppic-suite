#ifndef PROT_SPEC_THEO_PEAK_UTIL_HPP_
#define PROT_SPEC_THEO_PEAK_UTIL_HPP_

#include "spec/theo_peak.hpp"

namespace prot {

class TheoPeakUtil {
 public:
  static std::vector<double> getTheoMassVec (const TheoPeakPtrVec &theo_peak_list);
};

} /* namespace prot */

#endif 
