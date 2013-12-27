/*
 * comp_shift_low_mem.cpp
 *
 *  Created on: Dec 27, 2013
 *      Author: xunlikun
 */

#include <ptmsearch/comp_shift_low_mem.hpp>

namespace prot {

CompShiftLowMem::CompShiftLowMem(){
	for(unsigned int i=0;i< max_len_;i++){
		num_.push_back(0);
	}
}

} /* namespace prot */
