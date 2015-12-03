#include "base/logger.hpp"
#include "base/residue_util.hpp"

namespace prot {

static ResiduePtrVec convertStrToResiduePtrVec(const std::string &seq) {
  ResiduePtrVec residue_ptr_vec;
  for (size_t i = 0; i < seq.length(); i++) {
    AcidPtr acid_ptr = AcidBase::getAcidPtrByOneLetter(seq.substr(i, 1));
    PtmPtr ptm_ptr = PtmBase::getEmptyPtmPtr();
    ResiduePtr residue_ptr = ResidueBase::getBaseResiduePtr(acid_ptr, ptm_ptr);
    residue_ptr_vec.push_back(residue_ptr_vec);
  }
  return residue_ptr_vec;
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

}
