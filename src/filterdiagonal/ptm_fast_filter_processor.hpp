#ifndef PTM_FAST_FILTER_PROCESSOR_HPP_
#define PTM_FAST_FILTER_PROCESSOR_HPP_

#include "prsm/simple_prsm.hpp"
#include "filterdiagonal/ptm_fast_filter_mng.hpp"
#include "filterdiagonal/ptm_fast_filter_block.hpp"

namespace prot {

class PtmFastFilterProcessor {
public:
    PtmFastFilterProcessor(PtmFastFilterMngPtr mng_ptr);
    void process();

private:
    PtmFastFilterMngPtr mng_ptr_;
    PtmFastFilterBlockPtr filter_ptr_;

    void processBlock(int block, const std::string &sp_file_name,int n_spectra);
    SimplePrsmPtrVec2D readSimplePrsm(const std::string &sp_file_name);
    SimplePrsmPtrVec2D getBestMatch(const SimplePrsmPtrVec2D &match_ptrs,
                                    MsHeaderPtr header_ptr);
    void combineBlock();
};

typedef std::shared_ptr<PtmFastFilterProcessor> PtmFastFilterProcessorPtr;

} /* namespace prot */

#endif /* PTM_FAST_FILTER_PROCESSOR_HPP_ */
