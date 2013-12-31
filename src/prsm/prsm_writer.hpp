/*
 * prsm_writer.hpp
 *
 *  Created on: Dec 30, 2013
 *      Author: xunlikun
 */

#ifndef PRSM_WRITER_HPP_
#define PRSM_WRITER_HPP_

#include "prsm/prsm.hpp"

namespace prot {

class PrSMWriter {
public:
	int write(const char *prm_file_name);
	void addSimplePrSM(PrSMPtr matche);
	void addSimplePrSM(PrSMPtrVec matches);

private:
	PrSMPtrVec prsms_;
};

} /* namespace prot */

#endif /* PRSM_WRITER_HPP_ */
