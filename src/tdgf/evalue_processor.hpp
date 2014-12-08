#ifndef PROT_EVALUE_PROCESSOR_HPP_
#define PROT_EVALUE_PROCESSOR_HPP_

#include "htslib/faidx.h"

#include "base/proteoform.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_writer.hpp"
#include "tdgf/tdgf_mng.hpp"
#include "tdgf/comp_pvalue_array.hpp"
#include "tdgf/comp_pvalue_lookup_table.hpp"

namespace prot {

class EValueProcessor {
 public:
  EValueProcessor(TdgfMngPtr mng_ptr);
  ~EValueProcessor();

  void init();
  void process(bool is_separate);
  void processOneSpectrum(SpectrumSetPtr spec_set_ptr,
                          PrsmPtrVec &sele_prsm_ptrs,
                          bool is_separate,
                          PrsmWriter &writer);

  void compEvalues(SpectrumSetPtr spectrum_set_ptr, bool is_separate, 
                   PrsmPtrVec &sele_prsm_ptrs);

 private:
  TdgfMngPtr mng_ptr_;
  CompPValueArrayPtr comp_pvalue_ptr_;
  CompPValueLookupTablePtr comp_pvalue_table_ptr_;
  faidx_t *fai_;
  bool checkPrsms(const PrsmPtrVec &prsm_ptrs);

};

typedef std::shared_ptr<EValueProcessor> EValueProcessorPtr;

}

#endif
