/*
 * author  Xiaowen Liu
 * date    2013-11-01
 */

#ifndef PROT_BASE_RESIDUE_UTIL_HPP_
#define PROT_BASE_RESIDUE_UTIL_HPP_

#include <string>
#include <memory>
#include <map>

#include "base/prot_mod.hpp"
#include "base/residue.hpp"
#include "base/logger.hpp"

namespace prot {

class ResidueUtil {
 public:
  static ResiduePtrVec convertStrToResiduePtrVec(const std::string &seq);

  static ResiduePtrVec convertStrToResiduePtrVec(const std::string &seq, const ModPtrVec &fix_mod_ptr_vec);

  static ResiduePtrVec convertStrToResiduePtrVec(const StringPairVec &string_pair_vec);

  static ResiduePtrVec convertStrToResiduePtrVec(const StringPairVec &string_pair_vec,  
                                                 const ModPtrVec &fix_mod_ptr_vec);

  static int findResidue(const ResiduePtrVec &residue_list, ResiduePtr residue_ptr);

  static double compResiduePtrVecMass(const ResiduePtrVec &ptr_vec);
};

}
#endif
