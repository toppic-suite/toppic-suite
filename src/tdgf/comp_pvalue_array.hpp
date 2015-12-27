#ifndef PROT_COMP_PVALUE_ARRAY_HPP_
#define PROT_COMP_PVALUE_ARRAY_HPP_

#include "spec/spectrum_set.hpp"
#include "tdgf/tdgf_mng.hpp"
#include "tdgf/comp_prob_value.hpp"
#include "tdgf/count_test_num.hpp"

namespace prot {

class CompPValueArray {
 public:
  CompPValueArray(CountTestNumPtr test_num_ptr,
                  TdgfMngPtr mng_ptr);

  void compMultiExtremeValues(const PrmMsPtrVec &ms_six_ptr_vec, PrsmPtrVec &prsm_ptrs, 
                              double ppo, bool strict);

  void compSingleExtremeValue(const DeconvMsPtrVec &ms_ptr_vec, PrsmPtr prsm_ptr, double ppo);

  void process(SpectrumSetPtr spec_set_ptr, PrsmPtrVec &prsm_ptrs, double ppo, bool is_separate);

 private:
  TdgfMngPtr mng_ptr_;
  CompProbValuePtr comp_prob_ptr_;
  CountTestNumPtr test_num_ptr_;
  ResFreqPtrVec residue_ptrs_;
  ResFreqPtrVec pep_n_term_residue_ptrs_;
  ResFreqPtrVec prot_n_term_residue_ptrs_;

};

typedef std::shared_ptr<CompPValueArray> CompPValueArrayPtr;

}

#endif
