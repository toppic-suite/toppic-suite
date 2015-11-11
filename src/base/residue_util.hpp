/*
 * author  Xiaowen Liu
 * date    2013-11-01
 */

#ifndef PROT_BASE_RESIDUE_UTIL_HPP_
#define PROT_BASE_RESIDUE_UTIL_HPP_

#include <string>
#include <memory>
#include <map>

#include "base/residue.hpp"
#include "base/logger.hpp"

namespace prot {

class ResidueUtil {
 public:
  static ResiduePtr getResiduePtrByAcid(const ResiduePtrVec &residue_list,
                                        AcidPtr acid_ptr);

  static int findResidue(const ResiduePtrVec &residue_list, ResiduePtr residue_ptr);

  static ResiduePtrVec convertAcidToResiduePtrVec(const ResiduePtrVec &residue_list,
                                                  const AcidPtrVec &acid_list);

};

}
#endif
