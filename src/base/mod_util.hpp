#ifndef PROT_BASE_MOD_UTIL_HPP_
#define PROT_BASE_MOD_UTIL_HPP_

#include "base/mod.hpp"

namespace prot {

class ModUtil {
 public:
  static ModPtrVec readMod(const std::string &file_name);

  static ModPtrVec geneFixedModList(const std::string &str);

  static ResiduePtrVec geneResidueListWithMod(ResiduePtrVec residue_list,
                                              ModPtrVec fix_mod_list);
};

}
#endif
