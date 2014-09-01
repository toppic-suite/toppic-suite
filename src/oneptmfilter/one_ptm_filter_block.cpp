#include <iostream>
#include "oneptmfilter/one_ptm_filter_block.hpp"

namespace prot {

OnePtmFilterBlock::OnePtmFilterBlock(const ProteoformPtrVec &proteo_ptrs,
                                     OnePtmFilterMngPtr mng_ptr){
  mng_ptr_ = mng_ptr;
  proteo_ptrs_ = proteo_ptrs;
  LOG_DEBUG("start init blocks.")
  proteo_blocks_ = getProteoBlocks(proteo_ptrs, mng_ptr->db_block_size_);
  LOG_DEBUG("init blocks is done.")
}

void OnePtmFilterBlock::initBlock(int i) {
  filter_ptr_ = OnePtmFilterPtr(
      new OnePtmFilter(proteo_blocks_[i],mng_ptr_));
}

SimplePrsmPtrVec OnePtmFilterBlock::getBestMathBatch(
    SpectrumSetPtr spectrum_set_ptr){
  PrmMsPtr ms_ptr = spectrum_set_ptr->getMsTwoPtr();
  return filter_ptr_->getBestMatch(ms_ptr);
}

} /* namespace prot */
