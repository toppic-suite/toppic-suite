#ifndef PROT_BASE_TRUNC_UTIL_HPP_
#define PROT_BASE_TRUNC_UTIL_HPP_

#include "base/trunc.hpp"
#include "base/residue.hpp"

namespace prot {

class TruncUtil {
 public:
  static bool isValidTrunc(TruncPtr trunc_ptr, const ResiduePtrVec & res_ptr_vec);
};

}
#endif
