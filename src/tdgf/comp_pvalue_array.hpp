#ifndef PROT_COMP_PVALUE_ARRAY_HPP_
#define PROT_COMP_PVALUE_ARRAY_HPP_

#include "tdgf/tdgf_mng.hpp"
#include "tdgf/comp_prob_value.hpp"
#include "tdgf/count_test_num.hpp"

namespace prot {

class CompPValueArray {
 public:
  CompPValueArray(ProteoformPtrVec &raw_forms, 
                  ProteoformPtrVec &prot_mod_forms,
                  ResFreqPtrVec &residues,
                  TdgfMngPtr mng_ptr);

  ExtremeValuePtrVec compExtremeValue(PrmMsPtr ms_six, 
                                      PrSMPtrVec &prsms, bool strict);

 private:
  TdgfMngPtr mng_ptr_;
  CompProbValuePtr comp_prob_ptr_;
  CountTestNumPtr test_num_ptr_;
  ResFreqPtrVec pep_n_term_residues_;
  ResFreqPtrVec prot_n_term_residues_;
};


}

#endif
