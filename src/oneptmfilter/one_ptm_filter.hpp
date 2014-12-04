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
  void computeBestMatch(const PrmMsPtrVec &ms_ptr_vec);

  SimplePrsmPtrVec getCompMatchPtrs() {return comp_match_ptrs_;}
  SimplePrsmPtrVec getPrefMatchPtrs() {return pref_match_ptrs_;}
  SimplePrsmPtrVec getSuffMatchPtrs() {return suff_match_ptrs_;}
  SimplePrsmPtrVec getInternalMatchPtrs() {return internal_match_ptrs_;}

 private:
  OnePtmFilterMngPtr mng_ptr_;
  ProteoformPtrVec proteo_ptrs_;
  OnePtmCompShiftPtr index_ptr_;
  
  SimplePrsmPtrVec comp_match_ptrs_;
  SimplePrsmPtrVec pref_match_ptrs_;
  SimplePrsmPtrVec suff_match_ptrs_;
  SimplePrsmPtrVec internal_match_ptrs_;
};

typedef std::shared_ptr<OnePtmFilter> OnePtmFilterPtr;
} /* namespace prot */

#endif /* ONE_PTM_FILTER_H_ */
