#ifndef ONE_PTM_FILTER_PROCESSOR_HPP_
#define ONE_PTM_FILTER_PROCESSOR_HPP_

#include "prsm/simple_prsm.hpp"
#include "oneptmfilter/one_ptm_filter_mng.hpp"
#include "oneptmfilter/one_ptm_filter_block.hpp"

namespace prot {

class OnePtmFilterProcessor {
public:
    OnePtmFilterProcessor(OnePtmFilterMngPtr mng_ptr);
    void process();

private:
    OnePtmFilterMngPtr mng_ptr_;
    OnePtmFilterBlockPtr filter_ptr_;

    void processBlock(int block, const std::string &sp_file_name,int n_spectra);

};

typedef std::shared_ptr<PtmFastFilterProcessor> PtmFastFilterProcessorPtr;

} /* namespace prot */

#endif /* PTM_FAST_FILTER_PROCESSOR_HPP_ */
