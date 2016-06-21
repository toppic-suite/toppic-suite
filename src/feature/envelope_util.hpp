#ifndef PROT_FEATURE_ENVELOPE_UTIL_HPP_
#define PROT_FEATURE_ENVELOPE_UTIL_HPP_

#include <vector>

namespace prot {

class EnvelopeUtil {
 public:
  static int getMaxPos(std::vector<double> &values);
};

}

#endif
