#include "base/logger.hpp"
#include "base/residue_util.hpp"

namespace prot {

ResiduePtr ResidueUtil::getResiduePtrByAcid(const ResiduePtrVec &residue_ptr_vec,
                                            AcidPtr acid_ptr) {
  for (size_t i = 0; i < residue_ptr_vec.size(); i++) {
    if (residue_ptr_vec[i]->getAcidPtr() == acid_ptr) {
      return residue_ptr_vec[i];
    }
  }
  return ResiduePtr(nullptr);
}

int ResidueUtil::findResidue(const ResiduePtrVec &residue_list, 
                             ResiduePtr residue_ptr) {
  for (size_t i = 0; i < residue_list.size(); i++) {
    if (residue_list[i] == residue_ptr) {
      return i;
    }
  }
  return -1;
}


ResiduePtrVec ResidueUtil::convertAcidToResiduePtrVec(const ResiduePtrVec &residue_list,
                                                      const AcidPtrVec &acid_ptrs) {
  ResiduePtrVec result_seq;
  for (size_t i = 0; i < acid_ptrs.size(); i++) {
    ResiduePtr residue_ptr = getResiduePtrByAcid(residue_list, acid_ptrs[i]);
    result_seq.push_back(residue_ptr);
  }
  return result_seq;
}

}
