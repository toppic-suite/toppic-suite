/*
 * ptm_fast_filter_hi_mem.h
 *
 *  Created on: Dec 1, 2013
 *      Author: xunlikun
 */

#ifndef PTM_FAST_FILTER_HI_MEM_H_
#define PTM_FAST_FILTER_HI_MEM_H_

#include "ptm_fast_filter_mng.hpp"
#include "proteoform.hpp"

namespace prot {

class PtmFastFilterHiMem {
public:
private:
	PtmFastFilterMngPtr mng_;
	ProteoformPtrVec seqs_;

};

} /* namespace prot */

#endif /* PTM_FAST_FILTER_HI_MEM_H_ */
