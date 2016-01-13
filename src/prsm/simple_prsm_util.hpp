#ifndef PROT_PRSM_SIMPLE_PRSM_UTIL_HPP_
#define PROT_PRSM_SIMPLE_PRSM_UTIL_HPP_

#include "prsm/simple_prsm.hpp"

namespace prot {

class SimplePrsmUtil {
 public:
  static SimplePrsmPtrVec getUniqueMatches(SimplePrsmPtrVec &match_ptrs);
};

} /* namespace prot */

#endif /* SIMPLE_PRSM_HPP_ */
