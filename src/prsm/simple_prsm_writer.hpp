/*
 * simple_prsm_writer.hpp
 *
 *  Created on: Dec 9, 2013
 *      Author: xunlikun
 */

#ifndef PROT_SIMPLE_PRSM_WRITER_HPP_
#define PROT_SIMPLE_PRSM_WRITER_HPP_

#include "prsm/simple_prsm.hpp"

namespace prot {

class SimplePrSMWriter {
public:
	int write(const char *spectrum_file);
	void addSimplePrSM(SimplePrSMPtrVec matches);
private:
	SimplePrSMPtrVec matches_;
};

} /* namespace prot */

#endif /* SIMPLE_PRSM_WRITER_HPP_ */
