/*
 * ptm_fast_filter_hi_mem.h
 *
 *  Created on: Dec 1, 2013
 *      Author: xunlikun
 */

#ifndef PTM_FAST_FILTER_HI_MEM_H_
#define PTM_FAST_FILTER_HI_MEM_H_

#include "ptm_fast_filter_mng.hpp"
#include "base/proteoform.hpp"
#include "comp_shift_hi_mem.hpp"
#include "prsm/simple_prsm.hpp"
#include "spec/prm_peak.hpp"

namespace prot {

class PtmFastFilterHiMem {
public:
	PtmFastFilterHiMem(ProteoformPtrVec seqs,PtmFastFilterMngPtr mng);
	SimplePrSMPtrVec getBestMatch(PrmPeakMS ms);

private:
	PtmFastFilterMngPtr mng_;
	ProteoformPtrVec seqs_;
	CompShiftHiMemPtr index_;

	SimplePrSMPtrVec2D compute(PrmPeakMS ms);
	SimplePrSMPtrVec sort(SimplePrSMPtrVec2D matches);
};

typedef std::shared_ptr<PtmFastFilterHiMem> PtmFastFilterHiMemPtr;
} /* namespace prot */

#endif /* PTM_FAST_FILTER_HI_MEM_H_ */
