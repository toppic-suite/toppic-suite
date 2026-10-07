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

  // mod_residues: the PrSMs are variable PTM PrSMs and get the generating
  // function of the residues with their modified forms (var_residue_ptrs_)
  void compMultiExpectedValues(const PrmMsPtrVec& ms_six_ptr_vec,
                               PrsmPtrVec& prsm_ptrs, bool strict,
                               bool mod_residues = false);

  void compSingleExpectedValue(const DeconvMsPtrVec& ms_ptr_vec,
                               const PrsmPtr& prsm_ptr,
                               double n_term_label_mass);

  // is_separate false: the PrSMs of a spectrum share the generating
  // function of its PRM spectra, except those matched at +-1 isotope from
  // its precursor mass (getIsotopeOffset), which get the PRM spectra of
  // that mass, and, with the variable PTMs of TdgfMng::var_ptm_mod_ptr_vec_,
  // the variable PTM PrSMs, which get a generating function of their own
  void process(const SpectrumSetPtr& spec_set_ptr, PrsmPtrVec& prsm_ptrs,
               bool is_separate);

 private:
  // +1 or -1 for a PrSM without mass shifts whose adjusted precursor mass
  // differs from the spectrum's precursor mass by more than the precursor
  // error tolerance, 0 otherwise
  int getIsotopeOffset(const SpectrumSetPtr& spec_set_ptr,
                       const PrsmPtr& prsm_ptr);

  // a PrSM of the variable PTM search when the residues with their
  // modified forms are used: no unexpected shift, and variable PTMs or a
  // truncation N-terminal modification
  bool useModResidues(const PrsmPtr& prsm_ptr);

  TdgfMngPtr mng_ptr_;

  // the precursor error tolerance (ppm, or Dalton with toppic's
  // --prec-dalton-error): the candidate number's mass window, the isotope
  // offset and the zero and precursor PRM peaks of the generating function
  PeakTolerancePtr prec_tole_ptr_;

  CompProbValuePtr comp_prob_ptr_;

  CountTestNumPtr test_num_ptr_;

  ResFreqPtrVec residue_ptrs_;
  ResFreqPtrVec pep_n_term_residue_ptrs_;
  ResFreqPtrVec prot_n_term_residue_ptrs_;
  int residue_avg_len_ = 0;

  // the residues and their forms with the variable PTMs, for the variable
  // PTM PrSMs; empty when TdgfMng::var_ptm_mod_ptr_vec_ is
  ResFreqPtrVec var_residue_ptrs_;
  ResFreqPtrVec var_prot_n_term_residue_ptrs_;
  int var_residue_avg_len_ = 0;
};

using CompPValueArrayPtr = std::shared_ptr<CompPValueArray>;

}  // namespace toppic

#endif
