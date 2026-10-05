// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef TOPPIC_STAT_TDGF_COMP_PVALUE_ARRAY_HPP_
#define TOPPIC_STAT_TDGF_COMP_PVALUE_ARRAY_HPP_

#include "ms/spec/spectrum_set.hpp"
#include "stat/count/count_test_num.hpp"
#include "stat/tdgf/comp_prob_value.hpp"
#include "stat/tdgf/tdgf_mng.hpp"

namespace toppic {

class CompPValueArray {
 public:
  CompPValueArray(const CountTestNumPtr& test_num_ptr,
                  const TdgfMngPtr& mng_ptr);

  void compMultiExpectedValues(const PrmMsPtrVec& ms_six_ptr_vec,
                               PrsmPtrVec& prsm_ptrs, double ppo, bool strict);

  void compSingleExpectedValue(const DeconvMsPtrVec& ms_ptr_vec,
                               const PrsmPtr& prsm_ptr, double ppo,
                               double n_term_label_mass);

  // is_separate false: the PrSMs of a spectrum share the generating
  // function of its PRM spectra, except those matched at +-1 isotope from
  // its precursor mass (getIsotopeOffset), which get the PRM spectra of
  // that mass
  void process(const SpectrumSetPtr& spec_set_ptr, PrsmPtrVec& prsm_ptrs,
               double ppo, bool is_separate);

 private:
  // +1 or -1 for a PrSM without mass shifts whose adjusted precursor mass
  // differs from the spectrum's precursor mass by more than the precursor
  // error tolerance, 0 otherwise
  int getIsotopeOffset(const SpectrumSetPtr& spec_set_ptr,
                       const PrsmPtr& prsm_ptr, double ppo);

  TdgfMngPtr mng_ptr_;

  CompProbValuePtr comp_prob_ptr_;

  CountTestNumPtr test_num_ptr_;

  ResFreqPtrVec residue_ptrs_;
  ResFreqPtrVec pep_n_term_residue_ptrs_;
  ResFreqPtrVec prot_n_term_residue_ptrs_;
};

using CompPValueArrayPtr = std::shared_ptr<CompPValueArray>;

}  // namespace toppic

#endif
