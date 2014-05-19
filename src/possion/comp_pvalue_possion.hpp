#ifndef PROT_COMP_PVALUE_ARRAY_HPP_
#define PROT_COMP_PVALUE_ARRAY_HPP_

#include "tdgf/tdgf_mng.hpp"
#include "tdgf/count_test_num.hpp"

namespace prot {

class CompPValuePossion {
 public:
  CompPValuePossion(ProteoformPtrVec &raw_forms, 
                    ProteoformPtrVec &prot_mod_forms,
                    TdgfMngPtr mng_ptr);

  ExtremeValuePtrVec compExtremeValues(PrmMsPtr ms_six, 
                                       PrsmPtrVec &prsms, bool strict);

  void setPValue(DeconvMsPtr ms_ptr, PrsmPtr prsm_ptr);
  void setPValueArray(PrmMsPtr prm_ms_ptr, PrsmPtrVec prsms);

 private:
  TdgfMngPtr mng_ptr_;
  CountTestNumPtr test_num_ptr_;

  ExtremeValuePtr compExtremeValue(PrmMsPtr ms_ptr, PrsmPtr prsm_ptr);
};

typedef std::shared_ptr<CompPValuePossion> CompPValuePossionPtr;

}

#endif
