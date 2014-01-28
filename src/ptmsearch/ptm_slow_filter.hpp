/*
 * ptm_slow_filter.hpp
 *
 *  Created on: Dec 27, 2013
 *      Author: xunlikun
 */

#ifndef PROT_PTM_SLOW_FILTER_HPP_
#define PROT_PTM_SLOW_FILTER_HPP_

#include <memory>
#include "ptmsearch/comp_shift_low_mem.hpp"
#include "ptmsearch/ptm_slow_match.hpp"
#include "spec/spectrum_set.hpp"
#include "prsm/simple_prsm.hpp"
#include "ptmsearch/ptm_mng.hpp"

namespace prot {

class PtmSlowFilter {
public:
    PtmSlowFilter(
            SpectrumSetPtr spectrum_set,
            SimplePrSMPtrVec fast_Matches,
            CompShiftLowMemPtr comp_shift,
            PtmMngPtr mng);
    PtmSlowMatchPtrVec getBestMatch(int nshift,int type);
    PtmSlowMatchPtrVec getMatches(){return slow_matches_;};
private:
    PtmSlowMatchPtrVec slow_matches_;
};

typedef std::shared_ptr<PtmSlowFilter> PtmSlowFilterPtr;

} /* namespace prot */

#endif /* PTM_SLOW_FILTER_HPP_ */
