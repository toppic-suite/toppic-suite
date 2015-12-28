#ifndef PROT_ONE_PTM_SLOW_FILTER_HPP_
#define PROT_ONE_PTM_SLOW_FILTER_HPP_

#include <memory>

#include "spec/spectrum_set.hpp"
#include "prsm/simple_prsm.hpp"
#include "oneptmsearch/comp_shift_low_mem.hpp"
#include "oneptmsearch/one_ptm_search_mng.hpp"
#include "oneptmsearch/one_ptm_slow_match.hpp"

namespace prot {

class OnePtmSlowFilter {
 public:
  OnePtmSlowFilter(
      SpectrumSetPtr spectrum_set_ptr,
      SimplePrsmPtrVec simple_prsm_ptrs,
      CompShiftLowMemPtr comp_shift_ptr,
      SemiAlignTypePtr type_ptr, 
      OnePtmSearchMngPtr mng_ptr);
  PrsmPtrVec getPrsms() {return prsm_ptrs_;}

 private:
  PrsmPtrVec prsm_ptrs_;
};

typedef std::shared_ptr<OnePtmSlowFilter> OnePtmSlowFilterPtr;

} /* namespace prot */

#endif /* ONE_PTM_SLOW_FILTER_HPP_ */
