#include "prot_mod.hpp"

namespace prot {

ProtMod::ProtMod(std::string name, TruncPtr trunc_ptr, PtmPtr ptm_ptr,
                 AcidPtrVec valid_acid_ptrs) {
  name_ = name;
  trunc_ptr_ = trunc_ptr;
  ptm_ptr_ = ptm_ptr;
  valid_acid_ptrs_ = valid_acid_ptrs;
  prot_shift_ = trunc_ptr_->getShift() + ptm_ptr_->getMonoMass();
  pep_shift_ = ptm_ptr_->getMonoMass();
}

}
