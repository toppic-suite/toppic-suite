#include "base/logger.hpp"
#include "base/acid_base.hpp"
#include "base/ptm_base.hpp"
#include "base/residue_base.hpp"
#include "base/residue_util.hpp"

namespace prot {

ResiduePtrVec ResidueUtil::convertStrToResiduePtrVec(const std::string &seq) {
  ResiduePtrVec residue_ptr_vec;
  for (size_t i = 0; i < seq.length(); i++) {
    AcidPtr acid_ptr = AcidBase::getAcidPtrByOneLetter(seq.substr(i, 1));
    PtmPtr ptm_ptr = PtmBase::getEmptyPtmPtr();
    ResiduePtr residue_ptr = ResidueBase::getBaseResiduePtr(acid_ptr, ptm_ptr);
    residue_ptr_vec.push_back(residue_ptr);
  }
  return residue_ptr_vec;
}

void applyFixedMod(ResiduePtrVec &residue_ptrs, const ModPtrVec &fix_mod_ptr_vec) {
  for (size_t i = 0; i < residue_ptrs.size(); i++) {
    for (size_t j = 0; j < fix_mod_ptr_vec.size(); j++) {
      if (residue_ptrs[i] == fix_mod_ptr_vec[j]->getOriResiduePtr()) {
        residue_ptrs[i] = fix_mod_ptr_vec[j]->getModResiduePtr();
        break;
      }
    }
  }
}

ResiduePtrVec ResidueUtil::convertStrToResiduePtrVec(const std::string &seq, 
                                                     const ModPtrVec &fix_mod_ptr_vec) {
  ResiduePtrVec residue_ptrs = ResidueUtil::convertStrToResiduePtrVec(seq);
  applyFixedMod(residue_ptrs, fix_mod_ptr_vec);
  return residue_ptrs;
}

ResiduePtrVec ResidueUtil::convertStrToResiduePtrVec(const StringPairVec &string_pair_vec) {
  ResiduePtrVec residue_ptr_vec;
  for (size_t i = 0; i < string_pair_vec.size(); i++) {
    std::string acid_str = string_pair_vec[i].first;
    AcidPtr acid_ptr = AcidBase::getAcidPtrByOneLetter(acid_str);
    std::string ptm_str = string_pair_vec[i].second;
    PtmPtr ptm_ptr = PtmBase::getPtmPtrByAbbrName(ptm_str);
    ResiduePtr residue_ptr = ResidueBase::getBaseResiduePtr(acid_ptr, ptm_ptr);
    residue_ptr_vec.push_back(residue_ptr);
  }
  return residue_ptr_vec;
}

ResiduePtrVec ResidueUtil::convertStrToResiduePtrVec(const StringPairVec &string_pair_vec,  
                                                     const ModPtrVec &fix_mod_ptr_vec) {
  ResiduePtrVec residue_ptrs = ResidueUtil::convertStrToResiduePtrVec(string_pair_vec);
  applyFixedMod(residue_ptrs, fix_mod_ptr_vec);
  return residue_ptrs;
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

double ResidueUtil::compResiduePtrVecMass(const ResiduePtrVec &ptr_vec) {
  double mass = 0;
  for (size_t i = 0; i < ptr_vec.size(); i++) {
    mass += ptr_vec[i]->getMass();
  }
  return mass;
}

}
