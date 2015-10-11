#ifndef ZERO_PTM_FILTER_PROCESSOR_HPP_
#define ZERO_PTM_FILTER_PROCESSOR_HPP_

#include "base/db_block.hpp"
#include "prsm/simple_prsm.hpp"
#include "zeroptmfilter/zero_ptm_filter_mng.hpp"

namespace prot {

class ZeroPtmFilterProcessor {
public:
    ZeroPtmFilterProcessor(ZeroPtmFilterMngPtr mng_ptr);
    void process();

private:
    ZeroPtmFilterMngPtr mng_ptr_;

    void processBlock(DbBlockPtr block_ptr, int total_block_num);

};

typedef std::shared_ptr<ZeroPtmFilterProcessor> ZeroPtmFilterProcessorPtr;

} /* namespace prot */

#endif /* ZERO_PTM_FILTER_PROCESSOR_HPP_ */
