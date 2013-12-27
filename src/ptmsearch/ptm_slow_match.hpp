/*
 * ptm_slow_match.hpp
 *
 *  Created on: Dec 27, 2013
 *      Author: xunlikun
 */

#ifndef PTM_SLOW_MATCH_HPP_
#define PTM_SLOW_MATCH_HPP_

#include <memory>
#include <vector>
#include "prsm/prsm.hpp"

namespace prot {

class PtmSlowMatch {
public:
	double getScr(int shiftnum,int type);
	PrSMPtr geneResult(int shift_num,int type);
};

typedef std::shared_ptr<PtmSlowMatch> PtmSlowMatchPtr;
typedef std::vector<PtmSlowMatchPtr> PtmSlowMatchPtrVec;

} /* namespace prot */

#endif /* PTM_SLOW_MATCH_HPP_ */
