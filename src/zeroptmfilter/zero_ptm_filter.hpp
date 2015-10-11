#ifndef ZERO_PTM_FILTER_H_
#define ZERO_PTM_FILTER_H_

#include "base/proteoform.hpp"
#include "spec/extend_peak.hpp"
#include "prsm/simple_prsm.hpp"
#include "zeroptmfilter/zero_ptm_filter_mng.hpp"
#include "zeroptmfilter/zero_ptm_comp_shift.hpp"

namespace prot {

class ZeroPtmFilter {
 public:
  ZeroPtmFilter(const ProteoformPtrVec &proteo_ptrs,
                ZeroPtmFilterMngPtr mng_ptr);
  void computeBestMatch(const ExtendMsPtrVec &ms_ptr_vec);

  SimplePrsmPtrVec getCompMatchPtrs() {return comp_match_ptrs_;}
  SimplePrsmPtrVec getPrefMatchPtrs() {return pref_match_ptrs_;}
  SimplePrsmPtrVec getSuffMatchPtrs() {return suff_match_ptrs_;}
  SimplePrsmPtrVec getInternalMatchPtrs() {return internal_match_ptrs_;}

 private:
  ZeroPtmFilterMngPtr mng_ptr_;
  ProteoformPtrVec proteo_ptrs_;
  ZeroPtmCompShiftPtr index_ptr_;
  
  SimplePrsmPtrVec comp_match_ptrs_;
  SimplePrsmPtrVec pref_match_ptrs_;
  SimplePrsmPtrVec suff_match_ptrs_;
  SimplePrsmPtrVec internal_match_ptrs_;
};

typedef std::shared_ptr<ZeroPtmFilter> ZeroPtmFilterPtr;
} /* namespace prot */

#endif /* ZERO_PTM_FILTER_H_ */
