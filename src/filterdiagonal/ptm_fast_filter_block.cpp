#include <iostream>
#include "ptm_fast_filter_block.hpp"

namespace prot {

PtmFastFilterBlock::PtmFastFilterBlock(const ProteoformPtrVec &proteo_ptrs,
                                       PtmFastFilterMngPtr mng_ptr){
  mng_ptr_ = mng_ptr;
  proteo_ptrs_ = proteo_ptrs;
  LOG_DEBUG("start init blocks.")
  initProteoBlocks();
  LOG_DEBUG("init blocks is done.")
}

void PtmFastFilterBlock::initBlock(int i) {
  filter_ptr_ = PtmFastFilterHiMemPtr(
      new PtmFastFilterHiMem(proteo_blocks_[i],mng_ptr_));
}

void PtmFastFilterBlock::initProteoBlocks(){
  size_t start_idx = 0;
  size_t proteo_idx =0;
  int block_len=0;
  while(proteo_idx < proteo_ptrs_.size()){
    int proteo_len = proteo_ptrs_[proteo_idx]->getResSeqPtr()->getLen();
    if(block_len + proteo_len < mng_ptr_->db_block_size_) {
      block_len = block_len + proteo_len;
      proteo_idx++;
    }
    else{
      size_t end_idx = proteo_idx;
      ProteoformPtrVec proteo_in_block;
      for(size_t i = start_idx; i<=end_idx; i++){
        proteo_in_block.push_back(proteo_ptrs_[i]);
      }
      proteo_blocks_.push_back(proteo_in_block);
      start_idx = end_idx +1;
      proteo_idx = end_idx +1;
      block_len = 0;
    }
  }
  /* last block */
  if (start_idx < proteo_ptrs_.size()) {
    ProteoformPtrVec proteo_in_block;
    for(size_t i = start_idx; i < proteo_ptrs_.size(); i++){
      proteo_in_block.push_back(proteo_ptrs_[i]);
    }
    proteo_blocks_.push_back(proteo_in_block);
  }
}

SimplePrsmPtrVec PtmFastFilterBlock::getBestMathBatch(
    SpectrumSetPtr spectrum_set_ptr){
  SimplePrsmPtrVec result_ptrs;
  PrmMsPtr ms_ptr = spectrum_set_ptr->getMsTwoPtr();
  SimplePrsmPtrVec fast_match_list = filter_ptr_->getBestMatch(ms_ptr);
  for(size_t i=0; i<fast_match_list.size(); i++){
    if(i >= mng_ptr_->ptm_fast_filter_result_num_){
      break;
    }
    result_ptrs.push_back(fast_match_list[i]);
  }
  return result_ptrs;
}

} /* namespace prot */
