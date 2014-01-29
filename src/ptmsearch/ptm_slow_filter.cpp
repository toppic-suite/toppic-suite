/*
 * ptm_slow_filter.cpp
 *
 *  Created on: Dec 27, 2013
 *      Author: xunlikun
 */

#include "ptmsearch/ptm_slow_filter.hpp"

namespace prot {
PtmSlowFilter::PtmSlowFilter(
        SpectrumSetPtr spectrum_set,
        SimplePrSMPtrVec fast_Matches,
        CompShiftLowMemPtr comp_shift,
        PtmMngPtr mng){
    for(unsigned int i=0;i<fast_Matches.size();i++){
        ProteoformPtr seq = fast_Matches[i]->getSeq();
        slow_matches_.push_back(PtmSlowMatchPtr(
                new PtmSlowMatch(seq,spectrum_set,comp_shift,mng)));
    }
}
PtmSlowMatchPtrVec PtmSlowFilter::getBestMatch(int nshift,int type){
    PtmSlowMatchPtrVec matches;
    for(unsigned int i=0;i<slow_matches_.size();i++){
        slow_matches_[i]->setShift(nshift);
        slow_matches_[i]->setType(type);
        matches.push_back(slow_matches_[i]);
    }

    std::sort(matches.begin(),matches.end(),prot::psm_down);
    return matches;
}
} /* namespace prot */
