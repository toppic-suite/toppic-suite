#include "residue_util.hpp"
#include "xml_dom.hpp"
#include "xml_dom_document.hpp"

namespace proteomics {

ResiduePtr getResiduePtrByAcid(ResiduePtrVec residue_ptr_vec, 
                               AcidPtr acid_ptr) {
  for (unsigned int i = 0; i < residue_ptr_vec.size(); i++) {
    if (residue_ptr_vec[i]->getAcidPtr() == acid_ptr) {
      return residue_ptr_vec[i];
    }
  }
  return ResiduePtr(nullptr);
}


ResiduePtr getResiduePtrByAcidPtm(ResiduePtrVec residue_ptr_vec, 
                                  AcidPtr acid_ptr, PtmPtr ptm_ptr) {
  for (unsigned int i = 0; i < residue_ptr_vec.size(); i++) {
    if (residue_ptr_vec[i]->isSame(acid_ptr, ptm_ptr)) {
      return residue_ptr_vec[i];
    }
  }
  return ResiduePtr(nullptr);
}

}

