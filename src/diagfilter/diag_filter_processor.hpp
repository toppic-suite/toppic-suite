#ifndef PROT_DIAG_FILTER_PROCESSOR_HPP_
#define PROT_DIAG_FILTER_PROCESSOR_HPP_

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
    DiagFilterBlockPtr filter_ptr_;

    void processBlock(int block, const std::string &sp_file_name,int n_spectra);

};

typedef std::shared_ptr<DiagFilterProcessor> DiagFilterProcessorPtr;

} /* namespace prot */

#endif /* PROT_DIAG_FILTER_PROCESSOR_HPP_ */
