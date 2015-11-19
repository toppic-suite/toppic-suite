#include "base/logger.hpp"
#include "base/mod_base.hpp"
#include "base/prot_mod_base.hpp"
#include "base/prot_mod_util.hpp"
#include "base/trunc_util.hpp"

namespace prot {

bool ProtModUtil::allowMod(ProtModPtr prot_mod_ptr, const ResiduePtrVec &residues){
  if (prot_mod_ptr == ProtModBase::getProtModPtr_NONE()) {
    return true;
  }
  else {
    // check trunc
    if (!TruncUtil::isValidTrunc(prot_mod_ptr->getTruncPtr(), residues)) {
      return false;
    }
    ModPtr mod_ptr = prot_mod_ptr->getModPtr();
    if (mod_ptr != ModBase::getNoneModPtr()) {
      int mod_pos = prot_mod_ptr->getModPos();
      if (mod_pos < (int)residues.size()) { 
        return false;
      }
      if (residues[mod_pos] != mod_ptr->getOriResiduePtr()) {
        return false;
      }
    }
    return true;
  }
}

/*
bool ProtModUtil::contain_NME_ACETYLATION(const ProtModPtrVec &prot_mod_ptrs) {
  for (size_t i = 0; i < prot_mod_ptrs.size(); i++) {
    if (prot_mod_ptrs[i] == ProtModBase::getProtModPtr_NME_ACETYLATION()) {
      return true;
    }
  }
  return false;
}
*/

}
