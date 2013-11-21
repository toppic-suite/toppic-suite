/*
 * author  Xiaowen Liu
 * date    2013-11-1
 */
#include <stdlib.h>
#include <iostream>

#include "ptm.hpp"

namespace prot {

Ptm::Ptm(const std::string &abbr_name, 
            double mono_mass) {
    abbr_name_ = abbr_name;
    mono_mass_ = mono_mass;
}

bool Ptm::isEmpty() {
  if (mono_mass_ == 0.0) {
    return true;
  }
  else {
    return false;
  }
}
}

