#ifndef ZERO_PTM_FILTER_MASS_ZERO_PTM_FILTER_H_
#define ZERO_PTM_FILTER_MASS_ZERO_PTM_FILTER_H_

#include "base/proteoform.hpp"
#include "spec/extend_peak.hpp"
#include "spec/extend_ms.hpp"
#include "prsm/simple_prsm.hpp"
#include "zeroptmfilter/zero_ptm_filter_mng.hpp"
#include "zeroptmfilter/mass_match.hpp"

namespace prot {

class MassZeroPtmFilter {
 public:
  MassZeroPtmFilter(const ProteoformPtrVec &proteo_ptrs, ZeroPtmFilterMngPtr
                    mng_ptr);
  void computeBestMatch(const ExtendMsPtrVec &ms_ptr_vec);

  SimplePrsmPtrVec getCompMatchPtrs() {return comp_match_ptrs_;}
  SimplePrsmPtrVec getPrefMatchPtrs() {return pref_match_ptrs_;}
  SimplePrsmPtrVec getSuffMatchPtrs() {return suff_match_ptrs_;}
  SimplePrsmPtrVec getInternalMatchPtrs() {return internal_match_ptrs_;}

 private:
  ZeroPtmFilterMngPtr mng_ptr_;
  ProteoformPtrVec proteo_ptrs_;
  MassMatchPtr index_ptr_;
  
  SimplePrsmPtrVec comp_match_ptrs_;
  SimplePrsmPtrVec pref_match_ptrs_;
  SimplePrsmPtrVec suff_match_ptrs_;
  SimplePrsmPtrVec internal_match_ptrs_;
};

typedef std::shared_ptr<MassZeroPtmFilter> MassZeroPtmFilterPtr;
} /* namespace prot */

#endif /* ZERO_PTM_FILTER_H_ */
