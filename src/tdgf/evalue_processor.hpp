#ifndef PROT_EVALUE_PROCESSOR_HPP_
#define PROT_EVALUE_PROCESSOR_HPP_

#include "base/proteoform.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_writer.hpp"
#include "tdgf/tdgf_mng.hpp"
#include "tdgf/comp_pvalue_array.hpp"

namespace prot {

class EValueProcessor {
 public:
  EValueProcessor(TdgfMngPtr mng_ptr);
  void init();
  void process(bool is_separate);
  void processOneSpectrum(DeconvMsPtr ms_ptr, bool is_separate,
                          PrsmWriter &writer);
 private:
  TdgfMngPtr mng_ptr_;
  CompPValueArrayPtr comp_pvalue_ptr_;
  ProteoformPtrVec proteoforms_;
  PrsmPtrVec prsms_;
  bool checkPrsms(PrsmPtrVec &prsms);
};

}

#endif
