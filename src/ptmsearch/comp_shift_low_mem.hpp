/*
 * comp_shift_low_mem.hpp
 *
 *  Created on: Dec 27, 2013
 *      Author: xunlikun
 */

#ifndef COMP_SHIFT_LOW_MEM_HPP_
#define COMP_SHIFT_LOW_MEM_HPP_

#include <memory>
#include <vector>

namespace prot {

class CompShiftLowMem {
public:
	CompShiftLowMem();
	const static int max_len_ = 1000;


private:
	std::vector<short> num_;
};

typedef std::shared_ptr<CompShiftLowMem> CompShiftLowMemPtr;
} /* namespace prot */

#endif /* COMP_SHIFT_LOW_MEM_HPP_ */
