#ifndef PROT_DIAG_FILTER_H_
#define PROT_DIAG_FILTER_H_

#include "base/proteoform.hpp"
#include "spec/prm_ms.hpp"
#include "prsm/simple_prsm.hpp"
#include "zeroptmfilter/comp_shift.hpp"
#include "diagfilter/diag_filter_mng.hpp"

namespace prot {

class DiagFilter {
 public:
  DiagFilter(const ProteoformPtrVec &proteo_ptrs,
             DiagFilterMngPtr mng_ptr);
  SimplePrsmPtrVec getBestMatch(const PrmMsPtrVec &ms_ptr_vec);

 private:
  DiagFilterMngPtr mng_ptr_;
  ProteoformPtrVec proteo_ptrs_;
  CompShiftPtr index_ptr_;

  SimplePrsmPtrVec compute(const PrmMsPtrVec &ms_ptr_vec);
};

typedef std::shared_ptr<DiagFilter> DiagFilterPtr;
} /* namespace prot */

#endif /* PROT_DIAG_FILTER_H_ */
