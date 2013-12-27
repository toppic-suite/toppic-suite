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

namespace prot {

class PtmSearcher {
public:
	PtmSearcher(PtmMngPtr mng);
private:
	PtmMngPtr mng_;
//	CompShiftLowMemPtr
};

typedef std::shared_ptr<PtmSearcher> PtmSearcherPtr;
typedef std::vector<PtmSearcherPtr> PtmSearcherPtrVec;

} /* namespace prot */

#endif /* PTM_SEARCHER_HPP_ */
