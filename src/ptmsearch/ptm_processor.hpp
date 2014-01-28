/*
 * ptm_processor.hpp
 *
 *  Created on: Dec 20, 2013
 *      Author: xunlikun
 */

#ifndef PTM_PROCESSOR_HPP_
#define PTM_PROCESSOR_HPP_

#include "ptmsearch/ptm_mng.hpp"
#include "ptmsearch/ptm_searcher.hpp"
#include "base/proteoform.hpp"
#include "prsm/simple_prsm.hpp"

namespace prot {

class PtmProcessor {
public:
    PtmProcessor(PtmMngPtr mng);
    void process();
    void processDatabase(PtmSearcherPtr searcher);

    PtmMngPtr mng_;
    ProteoformPtrVec seqs_;
    SimplePrSMPtrVec simplePrsms_;
private:
    void init();
    void prsmFindSeq(SimplePrSMPtrVec simple_prsms,ProteoformPtrVec seqs);

};

typedef std::shared_ptr<PtmProcessor> PtmProcessorPtr;

} /* namespace prot */

#endif /* PTM_PROCESSOR_HPP_ */
