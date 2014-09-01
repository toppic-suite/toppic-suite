#ifndef PROT_DIAG_FILTER_H_
#define PROT_DIAG_FILTER_H_

#include "base/proteoform.hpp"
#include "spec/prm_peak.hpp"
#include "prsm/simple_prsm.hpp"
#include "diagfilter/diag_filter_mng.hpp"
#include "diagfilter/comp_shift.hpp"

namespace prot {

class DiagFilter {
 public:
  DiagFilter(const ProteoformPtrVec &proteo_ptrs,
             DiagFilterMngPtr mng_ptr);
  SimplePrsmPtrVec getBestMatch(PrmMsPtr ms_ptr);

 private:
  DiagFilterMngPtr mng_ptr_;
  ProteoformPtrVec proteo_ptrs_;
  CompShiftPtr index_ptr_;

  SimplePrsmPtrVec compute(PrmMsPtr ms_ptr);
};

typedef std::shared_ptr<DiagFilter> DiagFilterPtr;
} /* namespace prot */

#endif /* PROT_DIAG_FILTER_H_ */
