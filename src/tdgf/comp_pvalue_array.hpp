//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.


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

}  // namespace prot

#endif
