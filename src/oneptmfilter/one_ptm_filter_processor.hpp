#ifndef ONE_PTM_FILTER_PROCESSOR_HPP_
#define ONE_PTM_FILTER_PROCESSOR_HPP_

#include "base/db_block.hpp"
#include "prsm/simple_prsm.hpp"
#include "oneptmfilter/one_ptm_filter_mng.hpp"

namespace prot {

class OnePtmFilterProcessor {
public:
    OnePtmFilterProcessor(OnePtmFilterMngPtr mng_ptr);
    void process();

private:
    OnePtmFilterMngPtr mng_ptr_;

    void processBlock(DbBlockPtr block_ptr, int total_block_num);

};

typedef std::shared_ptr<OnePtmFilterProcessor> OnePtmFilterProcessorPtr;

} /* namespace prot */

#endif /* PTM_FAST_FILTER_PROCESSOR_HPP_ */
