/*
 * ptm_fast_filter_block.hpp
 *
 *  Created on: Dec 1, 2013
 *      Author: xunlikun
 */

#ifndef PTM_FAST_FILTER_BLOCK_HPP_
#define PTM_FAST_FILTER_BLOCK_HPP_

#include <memory>

#include "ptm_fast_filter_mng.hpp"
#include "proteoform.hpp"

namespace prot {

class PtmFastFilterBlock {
public:
	PtmFastFilterBlock();
private:
	PtmFastFilterMngPtr mng_;
	ProteoformPtrVec seqs_;
	std::vector<ProteoformPtrVec> seq_blocks;

};
typedef std::shared_ptr<PtmFastFilterBlock> PtmFastFilterBlockPtr;
} /* namespace prot */

#endif /* PTM_FAST_FILTER_BLOCK_HPP_ */
