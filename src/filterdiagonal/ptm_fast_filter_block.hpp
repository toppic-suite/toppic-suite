#ifndef PTM_FAST_FILTER_BLOCK_HPP_
#define PTM_FAST_FILTER_BLOCK_HPP_

#include <memory>

#include "base/proteoform.hpp"
#include "spec/spectrum_set.hpp"
#include "prsm/simple_prsm.hpp"
#include "filterdiagonal/ptm_fast_filter_mng.hpp"
#include "filterdiagonal/ptm_fast_filter_block.hpp"
#include "filterdiagonal/ptm_fast_filter_hi_mem.hpp"

namespace prot {

class PtmFastFilterBlock {
 public:
  PtmFastFilterBlock(ProteoformPtrVec seqs,PtmFastFilterMngPtr mng);
  int getBlockSize(){return seq_blocks_.size();}
  void initBlock(int i){
    filter_=NULL;
    filter_ = PtmFastFilterHiMemPtr(
        new PtmFastFilterHiMem(seq_blocks_[i],mng_));
  }
  SimplePrsmPtrVec getBestMathBatch(SpectrumSetPtr specttum_set);
 private:
  PtmFastFilterMngPtr mng_;
  ProteoformPtrVec seqs_;
  std::vector<ProteoformPtrVec> seq_blocks_;
  PtmFastFilterHiMemPtr filter_;

  void initSeqBlocks();
};

typedef std::shared_ptr<PtmFastFilterBlock> PtmFastFilterBlockPtr;
} /* namespace prot */

#endif /* PTM_FAST_FILTER_BLOCK_HPP_ */
