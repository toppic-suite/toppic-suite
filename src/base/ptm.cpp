/*
 * author  Xiaowen Liu
 * date    2013-11-1
 */
#include <stdlib.h>
#include <iostream>

#include "acid_vec.hpp"
#include "ptm.hpp"

namespace proteomics {

static std::string empty_ptm_name_ = "NON PTM";

Ptm::Ptm(const std::string &abbr_name, 
            const std::vector<AcidPtr> &valid_acid_ptrs, 
            double mono_mass) {
    abbr_name_ = abbr_name;
    valid_acid_ptrs_ = valid_acid_ptrs;
    mono_mass_ = mono_mass;
}

bool Ptm::isEmpty() {
  if (abbr_name_ == empty_ptm_name_) {
    return true;
  }
  else {
    return false;
  }
}

PtmPtr Ptm::getEmptyPtmPtr(std::vector<AcidPtr> &valid_acid_ptrs) {
  return PtmPtr(new Ptm(empty_ptm_name_, valid_acid_ptrs, 0));
}

}

