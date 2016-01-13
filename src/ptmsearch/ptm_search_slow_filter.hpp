#ifndef PROT_PTM_SEARCH_SLOW_FILTER_HPP_
#define PROT_PTM_SEARCH_SLOW_FILTER_HPP_

#include <memory>

#include "spec/spectrum_set.hpp"
#include "prsm/simple_prsm.hpp"
#include "oneptmsearch/ptm_search_mng.hpp"
#include "ptmsearch/ptm_slow_match.hpp"
#include "ptmsearch/comp_shift_low_mem.hpp"

namespace prot {

class PtmSearchSlowFilter {
 public:
  PtmSearchSlowFilter(
      SpectrumSetPtr spectrum_set_ptr,
      SimplePrsmPtrVec simple_prsm_ptrs,
      CompShiftLowMemPtr comp_shift_ptr,
      PtmSearchMngPtr mng_ptr);
  PrsmPtrVec getPrsms(int shift_num, AlignTypePtr type_ptr);

 private:
  PtmSlowMatchPtrVec complete_prefix_slow_match_ptrs_;
  PtmSlowMatchPtrVec suffix_internal_slow_match_ptrs_;
  PrsmPtrVec2D complete_prsm_2d_ptrs_;
  PrsmPtrVec2D prefix_prsm_2d_ptrs_;
  PrsmPtrVec2D suffix_prsm_2d_ptrs_;
  PrsmPtrVec2D internal_prsm_2d_ptrs_;
};

typedef std::shared_ptr<PtmSearchSlowFilter> PtmSearchSlowFilterPtr;

} /* namespace prot */

#endif /* PTM_SEARCH_SLOW_FILTER_HPP_ */
