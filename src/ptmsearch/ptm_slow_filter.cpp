/*
 * ptm_slow_filter.cpp
 *
 *  Created on: Dec 27, 2013
 *      Author: xunlikun
 */

#include <ptmsearch/ptm_slow_filter.hpp>

namespace prot {
PtmSlowFilter::PtmSlowFilter(SpectrumSetPtr spectrum_set,SimplePrSMPtrVec fast_Matches,CompShiftLowMemPtr comp_shift,PtmMngPtr mng){
	for(unsigned int i=0;i<fast_Matches.size();i++){
		ProteoformPtr seq = fast_Matches[i]->getSeq();
		slow_matches_.push_back(PtmSlowMatchPtr(new PtmSlowMatch()))
	}
}
} /* namespace prot */
