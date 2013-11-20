	
#include "residue.hpp"

namespace prot {

Residue::Residue(AcidPtr acid_ptr, PtmPtr ptm_ptr) {
  acid_ptr_ = acid_ptr;
  ptm_ptr_ = ptm_ptr;
  mass_ = acid_ptr->getMonoMass() + ptm_ptr->getMonoMass();
}

Residue::Residue(AcidPtrVec acid_ptr_vec, PtmPtrVec ptm_ptr_vec,
          std::string acid_name, std::string ptm_abbr_name) {
  acid_ptr_ = getAcidPtrByName(acid_ptr_vec, acid_name);
  ptm_ptr_ = getPtmPtrByAbbrName(ptm_ptr_vec, ptm_abbr_name);
  mass_ = acid_ptr_->getMonoMass() + ptm_ptr_->getMonoMass();
}

std::string Residue::toString(std::string delim_bgn, std::string delim_end) {
  if (ptm_ptr_->isEmpty()) {
    return acid_ptr_->getOneLetter();
  } else {
    return acid_ptr_->getOneLetter() + delim_bgn + ptm_ptr_->getAbbrName()
        + delim_end;
  }
}

}
