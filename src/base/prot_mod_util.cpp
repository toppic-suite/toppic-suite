#include "base/logger.hpp"
#include "base/prot_mod_base.hpp"
#include "base/prot_mod_util.hpp"

namespace prot {

bool ProtModUtil::allowMod(ProtModPtr prot_mod_ptr, const ResiduePtrVec &residues){
  if (prot_mod_ptr == ProtModBase::getProtModPtr_NONE()) {
    return true;
  }
  else if (prot_mod_ptr == ProtModBase::getProtModPtr_NME()) {
    if(residues.size()>=2 && residues[0]->getAcidPtr()->getOneLetter() == "M"){
      return true;
    }
    return false;
  }
  else if (prot_mod_ptr == ProtModBase::getProtModPtr_NME_ACETYLATION()){
    if(residues.size()>=2 && residues[0]->getAcidPtr()->getOneLetter() == "M"){
      return true;
    }
    return false;
  }
  return false;
}

bool ProtModUtil::contain_NME_ACETYLATION(const ProtModPtrVec &prot_mod_ptrs) {
  for (size_t i = 0; i < prot_mod_ptrs.size(); i++) {
    if (prot_mod_ptrs[i] == ProtModBase::getProtModPtr_NME_ACETYLATION()) {
      return true;
    }
  }
  return false;
}

}
