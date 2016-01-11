#ifndef PROT_SPEC_EXTEND_MS_HPP_
#define PROT_SPEC_EXTEND_MS_HPP_

#include <vector>
#include <memory>

#include "spec/extend_peak.hpp"
#include "spec/ms.hpp"

namespace prot {

typedef std::shared_ptr<Ms<ExtendPeakPtr>> ExtendMsPtr;
typedef std::vector<ExtendMsPtr> ExtendMsPtrVec;

class ExtendMs {
 public:
  static std::vector<double> getExtendMassVec(ExtendMsPtr extend_ms_ptr);

  static std::vector<std::pair<int, int>> getExtendIntMassErrorList(
      const ExtendMsPtrVec &ext_ms_ptr_vec, bool pref, double scale);

};

} /* namespace prot */

#endif 
