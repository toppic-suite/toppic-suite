#ifndef PROT_DIAG_FILTER_PROCESSOR_HPP_
#define PROT_DIAG_FILTER_PROCESSOR_HPP_

#include "base/db_block.hpp"
#include "prsm/simple_prsm.hpp"
#include "diagfilter/diag_filter_mng.hpp"
#include "diagfilter/diag_filter_block.hpp"

namespace prot {

class DiagFilterProcessor {
public:
    DiagFilterProcessor(DiagFilterMngPtr mng_ptr);
    void process();

private:
    DiagFilterMngPtr mng_ptr_;

    void processBlock(DbBlockPtr block_ptr, int total_block_num);

};

typedef std::shared_ptr<DiagFilterProcessor> DiagFilterProcessorPtr;

} /* namespace prot */

#endif /* PROT_DIAG_FILTER_PROCESSOR_HPP_ */
