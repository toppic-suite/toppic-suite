#ifndef ONE_PTM_FILTER_H_
#define ONE_PTM_FILTER_H_

#include "base/proteoform.hpp"
#include "spec/prm_peak.hpp"
#include "prsm/simple_prsm.hpp"
#include "oneptmfilter/one_ptm_filter_mng.hpp"
#include "oneptmfilter/one_ptm_comp_shift.hpp"

namespace prot {

class OnePtmFilter {
 public:
  OnePtmFilter(const ProteoformPtrVec &proteo_ptrs,
               OnePtmFilterMngPtr mng_ptr);
  SimplePrsmPtrVec getBestMatch(PrmMsPtr ms_ptr);

 private:
  OnePtmFilterMngPtr mng_ptr_;
  ProteoformPtrVec proteo_ptrs_;
  OnePtmCompShiftPtr index_ptr_;

  SimplePrsmPtrVec compute(PrmMsPtr ms_ptr);
};

typedef std::shared_ptr<OnePtmFilter> OnePtmFilterPtr;
} /* namespace prot */

#endif /* ONE_PTM_FILTER_H_ */
