/*
 * author  Xiaowen Liu
 * date    2013-11-1
 */
#include <stdlib.h>
#include <iostream>

#include "ptm.hpp"

namespace proteomics {

Ptm::Ptm(const std::string &abbr_name, 
            const std::vector<Acid> &valid_acids, 
            double mono_mass, 
            bool is_empty) {
    abbr_name_ = abbr_name;
    valid_acids_ = valid_acids;
    mono_mass_ = mono_mass;
    is_empty_ = is_empty;
}


}

