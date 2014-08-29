#ifndef PROT_COMP_PVALUE_ARRAY_HPP_
#define PROT_COMP_PVALUE_ARRAY_HPP_

#include "spec/spectrum_set.hpp"
#include "tdgf/tdgf_mng.hpp"
#include "tdgf/comp_prob_value.hpp"
#include "tdgf/count_test_num.hpp"

namespace prot {

class CompPValueArray {
 public:
  CompPValueArray(const ProteoformPtrVec &raw_proteo_ptrs, 
                  const ProteoformPtrVec &mod_proteo_ptrs,
                  const ResFreqPtrVec &residue_ptrs,
                  TdgfMngPtr mng_ptr);

  void compMultiExtremeValues(PrmMsPtr ms_six_ptr, PrsmPtrVec &prsm_ptrs, 
                              bool strict);

  void compSingleExtremeValue(DeconvMsPtr deconv_ms_ptr, PrsmPtr prsm_ptr);

  void process(SpectrumSetPtr spec_set_ptr, bool is_separate, PrsmPtrVec &prsm_ptrs);

 private:
  TdgfMngPtr mng_ptr_;
  CompProbValuePtr comp_prob_ptr_;
  CountTestNumPtr test_num_ptr_;
  ResFreqPtrVec pep_n_term_residue_ptrs_;
  ResFreqPtrVec prot_n_term_residue_ptrs_;

};

typedef std::shared_ptr<CompPValueArray> CompPValueArrayPtr;

}

#endif
