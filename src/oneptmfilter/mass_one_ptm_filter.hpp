#ifndef PROT_MASS_ONE_PTM_FILTER_H_
#define PROT_MASS_ONE_PTM_FILTER_H_

#include "base/proteoform.hpp"
#include "spec/prm_ms.hpp"
#include "prsm/simple_prsm.hpp"
#include "zeroptmfilter/mass_match.hpp"
#include "oneptmfilter/one_ptm_filter_mng.hpp"

namespace prot {

class MassOnePtmFilter {
 public:
  MassOnePtmFilter(const ProteoformPtrVec &proteo_ptrs,
               OnePtmFilterMngPtr mng_ptr);
  void computeBestMatch(const PrmMsPtrVec &prm_ms_ptr_vec,
                        const PrmMsPtrVec &srm_ms_ptr_vec);

  SimplePrsmPtrVec getCompMatchPtrs() {return comp_match_ptrs_;}
  SimplePrsmPtrVec getPrefMatchPtrs() {return pref_match_ptrs_;}
  SimplePrsmPtrVec getSuffMatchPtrs() {return suff_match_ptrs_;}
  SimplePrsmPtrVec getInternalMatchPtrs() {return internal_match_ptrs_;}

 private:
  OnePtmFilterMngPtr mng_ptr_;
  ProteoformPtrVec proteo_ptrs_;

  MassMatchPtr diag_index_ptr_;
  MassMatchPtr rev_diag_index_ptr_;
  MassMatchPtr term_index_ptr_;
  MassMatchPtr rev_term_index_ptr_;
  
  SimplePrsmPtrVec comp_match_ptrs_;
  SimplePrsmPtrVec pref_match_ptrs_;
  SimplePrsmPtrVec suff_match_ptrs_;
  SimplePrsmPtrVec internal_match_ptrs_;
};

typedef std::shared_ptr<MassOnePtmFilter> MassOnePtmFilterPtr;
} /* namespace prot */

#endif /* ONE_PTM_FILTER_H_ */
