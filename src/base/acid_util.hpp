#ifndef PROT_BASE_ACID_UTIL_HPP_
#define PROT_BASE_ACID_UTIL_HPP_

#include "base/acid.hpp"
#include "base/acid_base.hpp"

namespace prot {

class AcidUtil {
 public:
  /**
   * Converts a protein sequence (with one letter representation of amino
   * acids) to an amino acid array.
   */
  static AcidPtrVec convertStrToAcidPtrVec(const std::string &seq);

};

}
#endif
