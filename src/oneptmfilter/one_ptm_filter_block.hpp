#ifndef ONE_PTM_FILTER_BLOCK_HPP_
#define ONE_PTM_FILTER_BLOCK_HPP_

#include <memory>

#include "base/proteoform.hpp"
#include "spec/spectrum_set.hpp"
#include "prsm/simple_prsm.hpp"
#include "oneptmfilter/one_ptm_filter_mng.hpp"
#include "oneptmfilter/one_ptm_filter_block.hpp"
#include "oneptmfilter/one_ptm_filter.hpp"

namespace prot {

class OnePtmFilterBlock {
 public:
  OnePtmFilterBlock(const ProteoformPtrVec &proteo_ptrs, 
                    OnePtmFilterMngPtr mng_ptr);

  int getBlockSize(){return proteo_blocks_.size();}

  void initBlock(int i);

  SimplePrsmPtrVec getBestMathBatch(SpectrumSetPtr specttum_set_ptr);

 private:
  OnePtmFilterMngPtr mng_ptr_;
  ProteoformPtrVec proteo_ptrs_;
  ProteoformPtrVec2D proteo_blocks_;
  OnePtmFilterPtr filter_ptr_;

};

typedef std::shared_ptr<OnePtmFilterBlock> OnePtmFilterBlockPtr;
} /* namespace prot */

#endif /* ONE_PTM_FILTER_BLOCK_HPP_ */
