#ifndef PROT_EVALUE_PROCESSOR_HPP_
#define PROT_EVALUE_PROCESSOR_HPP_

#include "base/proteoform.hpp"
#include "prsm/prsm.hpp"
#include "tdgf/tdgf_mng.hpp"
#include "tdgf/comp_pvalue_array.hpp"

namespace prot {

class EValueProcessor {
 public:
  EValueProcessor(TdgfMngPtr mng_ptr);
  void init();

 private:
    TdgfMngPtr mng_ptr_;
    CompPValueArrayPtr comp_pvalue_ptr_;
    ProteoformPtrVec proteoforms_;
    PrSMPtrVec prsms_;
};

}

#endif
