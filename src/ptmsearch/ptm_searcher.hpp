/*
 * ptm_searcher.hpp
 *
 *  Created on: Dec 23, 2013
 *      Author: xunlikun
 */

#ifndef PTM_SEARCHER_HPP_
#define PTM_SEARCHER_HPP_

#include <memory>
#include <vector>
#include "ptmsearch/ptm_mng.hpp"
#include "spec/spectrum_set.hpp"
#include "prsm/simple_prsm.hpp"
#include "prsm/prsm.hpp"
#include "ptmsearch/ptm_slow_filter.hpp"
#include "ptmsearch/comp_shift_low_mem.hpp"

namespace prot {

class PtmSearcher {
public:
	PtmSearcher(PtmMngPtr mng);
	void search(SpectrumSetPtr spectrum_set,SimplePrSMPtrVec matches,PrSMPtrVec3D &prsms);
	void search(SpectrumSetPtr spectrum_set,SimplePrSMPtrVec matches,PrSMPtrVec3D &prsms,PtmSlowFilterPtr slow_filter);
private:
	PtmMngPtr mng_;
	CompShiftLowMemPtr comp_shift_;
};

typedef std::shared_ptr<PtmSearcher> PtmSearcherPtr;
typedef std::vector<PtmSearcherPtr> PtmSearcherPtrVec;

} /* namespace prot */

#endif /* PTM_SEARCHER_HPP_ */
