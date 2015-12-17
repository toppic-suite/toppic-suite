#ifndef PROT_PTM_SLOW_FILTER_HPP_
#define PROT_PTM_SLOW_FILTER_HPP_

#include <memory>

#include "spec/spectrum_set.hpp"
#include "prsm/simple_prsm.hpp"
#include "ptmsearch/ptm_mng.hpp"
#include "ptmsearch/comp_shift_low_mem.hpp"
#include "ptmsearch/ptm_slow_match.hpp"

namespace prot {

class PtmSlowFilter {
 public:
  PtmSlowFilter(
      SpectrumSetPtr spectrum_set_ptr,
      SimplePrsmPtrVec simple_prsm_ptrs,
      CompShiftLowMemPtr comp_shift_ptr,
      PtmMngPtr mng_ptr);
  PrsmPtrVec getPrsms(int shift_num, AlignTypePtr type_ptr);

 private:
  PtmSlowMatchPtrVec complete_prefix_slow_match_ptrs_;
  PtmSlowMatchPtrVec suffix_internal_slow_match_ptrs_;
  PrsmPtrVec2D complete_prsm_2d_ptrs_;
  PrsmPtrVec2D prefix_prsm_2d_ptrs_;
  PrsmPtrVec2D suffix_prsm_2d_ptrs_;
  PrsmPtrVec2D internal_prsm_2d_ptrs_;
};

typedef std::shared_ptr<PtmSlowFilter> PtmSlowFilterPtr;

} /* namespace prot */

#endif /* PTM_SLOW_FILTER_HPP_ */
