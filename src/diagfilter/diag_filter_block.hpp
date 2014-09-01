#ifndef PROT_DIAG_FILTER_BLOCK_HPP_
#define PROT_DIAG_FILTER_BLOCK_HPP_

#include <memory>

#include "base/proteoform.hpp"
#include "spec/spectrum_set.hpp"
#include "prsm/simple_prsm.hpp"
#include "diagfilter/diag_filter_mng.hpp"
#include "diagfilter/diag_filter_block.hpp"
#include "diagfilter/diag_filter.hpp"

namespace prot {

class DiagFilterBlock {
 public:
  DiagFilterBlock(const ProteoformPtrVec &proteo_ptrs, 
                  DiagFilterMngPtr mng_ptr);

  int getBlockSize(){return proteo_blocks_.size();}

  void initBlock(int i);

  SimplePrsmPtrVec getBestMathBatch(SpectrumSetPtr specttum_set_ptr);

 private:
  DiagFilterMngPtr mng_ptr_;
  ProteoformPtrVec proteo_ptrs_;
  ProteoformPtrVec2D proteo_blocks_;
  DiagFilterPtr filter_ptr_;
};

typedef std::shared_ptr<DiagFilterBlock> DiagFilterBlockPtr;
} /* namespace prot */

#endif /* PROT_DIAG_FILTER_BLOCK_HPP_ */
