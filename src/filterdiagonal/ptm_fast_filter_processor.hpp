/*
 * ptm_fast_filter_processor.hpp
 *
 *  Created on: Dec 1, 2013
 *      Author: xunlikun
 */

#ifndef PTM_FAST_FILTER_PROCESSOR_HPP_
#define PTM_FAST_FILTER_PROCESSOR_HPP_

#include "ptm_fast_filter_mng.hpp"
#include "ptm_fast_filter_block.hpp"
#include "prsm/simple_prsm.hpp"

namespace prot {

class PtmFastFilterProcessor {
public:
    PtmFastFilterProcessor(PtmFastFilterMngPtr mng);
    void process();

private:
    PtmFastFilterMngPtr mng_ptr_;
    PtmFastFilterBlockPtr filter_ptr_;

    void processBlock(int block,std::string sp_file_name,int n_spectra);
    void combineBlock(std::string sp_file_name);
    SimplePrSMPtrVec2D readSimplePrSM(std::string sp_file_name);
    SimplePrSMPtrVec2D getBestMatch(SimplePrSMPtrVec2D matches,
                                  MsHeaderPtr header);
};

typedef std::shared_ptr<PtmFastFilterProcessor> PtmFastFilterProcessorPtr;

} /* namespace prot */

#endif /* PTM_FAST_FILTER_PROCESSOR_HPP_ */
