#ifndef PTM_FAST_FILTER_HI_MEM_H_
#define PTM_FAST_FILTER_HI_MEM_H_

#include "base/proteoform.hpp"
#include "spec/prm_peak.hpp"
#include "prsm/simple_prsm.hpp"
#include "filterdiagonal/ptm_fast_filter_mng.hpp"
#include "filterdiagonal/comp_shift_hi_mem.hpp"

namespace prot {

class PtmFastFilterHiMem {
 public:
  PtmFastFilterHiMem(const ProteoformPtrVec &proteo_ptrs,
                     PtmFastFilterMngPtr mng_ptr);
  SimplePrsmPtrVec getBestMatch(PrmMsPtr ms_ptr);

 private:
  PtmFastFilterMngPtr mng_ptr_;
  ProteoformPtrVec proteo_ptrs_;
  CompShiftHiMemPtr index_ptr_;

  SimplePrsmPtrVec2D compute(PrmMsPtr ms_ptr);
  SimplePrsmPtrVec sort(const SimplePrsmPtrVec2D &match_ptrs);
};

typedef std::shared_ptr<PtmFastFilterHiMem> PtmFastFilterHiMemPtr;
} /* namespace prot */

#endif /* PTM_FAST_FILTER_HI_MEM_H_ */
