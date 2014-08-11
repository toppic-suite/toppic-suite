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
  PtmFastFilterBlock(const ProteoformPtrVec &proteo_ptrs, 
                     PtmFastFilterMngPtr mng_ptr);

  int getBlockSize(){return proteo_blocks_.size();}

  void initBlock(int i);

  SimplePrsmPtrVec getBestMathBatch(SpectrumSetPtr specttum_set_ptr);

 private:
  PtmFastFilterMngPtr mng_ptr_;
  ProteoformPtrVec proteo_ptrs_;
  std::vector<ProteoformPtrVec> proteo_blocks_;
  PtmFastFilterHiMemPtr filter_ptr_;

  void initProteoBlocks();
};

typedef std::shared_ptr<PtmFastFilterBlock> PtmFastFilterBlockPtr;
} /* namespace prot */

#endif /* PTM_FAST_FILTER_BLOCK_HPP_ */
