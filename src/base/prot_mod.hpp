#ifndef PROT_PROT_MOD_H_
#define PROT_PROT_MOD_H_

#include "ptm.hpp"
#include "trunc.hpp"

namespace prot {

class ProtMod {
 public:
  ProtMod(std::string name, TruncPtr trunc_ptr, PtmPtr ptm_ptr,
          AcidPtrVec valid_acid_ptrs);

  std::string getName() {return name_;};

  TruncPtr getTruncPtr() {return trunc_ptr_;}

  double getProtShift() {return prot_shift_;}

  double getPepShift() {return pep_shift_;}

 private:
  std::string name_;
  TruncPtr trunc_ptr_;
  PtmPtr ptm_ptr_;
  AcidPtrVec valid_acid_ptrs_;
  double prot_shift_;
  double pep_shift_;
};

}
#endif
