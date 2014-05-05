/*
 * ptm_slow_filter.hpp
 *
 *  Created on: Dec 27, 2013
 *      Author: xunlikun
 */

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
      SpectrumSetPtr spectrum_set,
      SimplePrsmPtrVec fast_Matches,
      CompShiftLowMemPtr comp_shift,
      PtmMngPtr mng);
  PrsmPtrVec getPrsms(int nshift, SemiAlignTypePtr type);
 private:
  PtmSlowMatchPtrVec complete_prefix_slow_matches_;
  PtmSlowMatchPtrVec suffix_internal_slow_matches_;
  PrsmPtrVec2D complete_prsms_;
  PrsmPtrVec2D prefix_prsms_;
  PrsmPtrVec2D suffix_prsms_;
  PrsmPtrVec2D internal_prsms_;
};

typedef std::shared_ptr<PtmSlowFilter> PtmSlowFilterPtr;

} /* namespace prot */

#endif /* PTM_SLOW_FILTER_HPP_ */
