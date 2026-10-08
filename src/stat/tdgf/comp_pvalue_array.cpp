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

#include "stat/tdgf/comp_pvalue_array.hpp"

#include <cmath>
#include <map>

#include "common/base/proteoform_type.hpp"
#include "common/base/mass_constant.hpp"
#include "common/util/logger.hpp"
#include "ms/factory/prm_ms_factory.hpp"
#include "ms/spec/deconv_ms_util.hpp"
#include "stat/count/count_util.hpp"
#include "stat/tdgf/comp_prob_value_array.hpp"

namespace toppic {

CompPValueArray::CompPValueArray(const CountTestNumPtr& test_num_ptr,
                                 const TdgfMngPtr& mng_ptr) {
  mng_ptr_ = mng_ptr;
  prec_tole_ptr_ =
      mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getPrecTolerancePtr();
  test_num_ptr_ = test_num_ptr;
  residue_ptrs_ = test_num_ptr->getResFreqPtrVec();
  pep_n_term_residue_ptrs_ = residue_ptrs_;
  prot_n_term_residue_ptrs_ = test_num_ptr->getNTermResFreqPtrVec();
  residue_avg_len_ = static_cast<int>(test_num_ptr->getResidueAvgLen());
  comp_prob_ptr_ = std::make_shared<CompProbValue>(
      mng_ptr_->convert_ratio_, residue_ptrs_, test_num_ptr->getResidueAvgLen(),
      mng_ptr_->unexpected_shift_num_ + 1, mng_ptr_->max_table_height_,
      mng_ptr_->max_prec_mass_,
      mng_ptr_->prsm_para_ptr_->useEvalueStrictComplement());
  // the variable PTM PrSMs' residues: the residues and their modified forms,
  // which share the residue's frequency (toppic_test experiment 129)
  if (!mng_ptr_->var_ptm_mod_ptr_vec_.empty()) {
    var_residue_ptrs_ = count_util::geneModResFreqs(
        residue_ptrs_, mng_ptr_->var_ptm_mod_ptr_vec_);
    var_prot_n_term_residue_ptrs_ = count_util::geneModResFreqs(
        prot_n_term_residue_ptrs_, mng_ptr_->var_ptm_mod_ptr_vec_);
    var_residue_avg_len_ = count_util::computeAvgLength(
        var_residue_ptrs_, mng_ptr_->convert_ratio_);
  }
}

bool CompPValueArray::useModResidues(const PrsmPtr& prsm_ptr) {
  if (var_residue_ptrs_.empty()) {
    return false;
  }
  ProteoformPtr proteo_ptr = prsm_ptr->getProteoformPtr();
  return proteo_ptr->getAlterNum(AlterType::UNEXPECTED) == 0 &&
         (proteo_ptr->getAlterNum(AlterType::VARIABLE) > 0 ||
          proteo_ptr->hasTruncNTermMod());
}

// set alignment
void CompPValueArray::compMultiExpectedValues(const PrmMsPtrVec& ms_six_ptr_vec,
                                              PrsmPtrVec& prsm_ptrs,
                                              bool strict, bool mod_residues) {
  PrmPeakPtrVec2D prm_ptr_2d;
  for (size_t i = 0; i < ms_six_ptr_vec.size(); i++) {
    PrmPeakPtrVec prm_peak_ptrs = ms_six_ptr_vec[i]->getPeakPtrVec();
    if (prm_peak_ptrs.size() > 0) {
      prm_ptr_2d.push_back(prm_peak_ptrs);
    }
  }
  double prob_prec_mass;
  // use proteoform mass when the strict mode is used, otherwise, the precursor
  // mass
  if (strict) {
    prob_prec_mass = prsm_ptrs[0]->getProteoformPtr()->getMass();
  } else {
    prob_prec_mass =
        ms_six_ptr_vec[0]->getMsHeaderPtr()->getFirstPrecMonoMass();
  }

  // the protein probabilities are used only by COMPLETE and PREFIX PrSMs
  // and the peptide ones only by SUFFIX and INTERNAL PrSMs, so a spectrum
  // whose PrSMs are all of one kind needs only one of the two generating
  // function computations
  bool prot_needed = false;
  bool pep_needed = false;
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    ProteoformTypePtr type_ptr =
        prsm_ptrs[i]->getProteoformPtr()->getProteoformType();
    if (type_ptr == ProteoformType::COMPLETE ||
        type_ptr == ProteoformType::PREFIX) {
      prot_needed = true;
    } else {
      pep_needed = true;
    }
  }
  // the residues of the generating function: the variable PTM PrSMs' with
  // mod_residues, otherwise the residues without variable PTMs; a table
  // serves both, so it is told which before every computation
  if (!var_residue_ptrs_.empty()) {
    comp_prob_ptr_->setResidues(
        mod_residues ? var_residue_ptrs_ : residue_ptrs_,
        mod_residues ? var_residue_avg_len_ : residue_avg_len_);
  }
  const ResFreqPtrVec& prot_n_term_residue_ptrs =
      mod_residues ? var_prot_n_term_residue_ptrs_ : prot_n_term_residue_ptrs_;
  const ResFreqPtrVec& pep_n_term_residue_ptrs =
      mod_residues ? var_residue_ptrs_ : pep_n_term_residue_ptrs_;
  std::vector<double> prot_probs;
  if (prot_needed) {
    comp_prob_value_array::compProbArray(
        comp_prob_ptr_, prot_n_term_residue_ptrs, prm_ptr_2d, prsm_ptrs,
        strict, prob_prec_mass, prec_tole_ptr_, prot_probs);
  }
  std::vector<double> pep_probs;
  if (pep_needed) {
    comp_prob_value_array::compProbArray(
        comp_prob_ptr_, pep_n_term_residue_ptrs, prm_ptr_2d, prsm_ptrs,
        strict, prob_prec_mass, prec_tole_ptr_, pep_probs);
  }

  double tolerance = prec_tole_ptr_->compPrecErrorTole(
      ms_six_ptr_vec[0]->getMsHeaderPtr()->getFirstPrecMonoMass());
  // the precursor mass of the PRM spectra (label-free): the original
  // precursor's, or the adjusted one for the PrSMs of a +-1 isotope offset
  // (process)
  double spec_prec_mass =
      ms_six_ptr_vec[0]->getMsHeaderPtr()->getFirstPrecMonoMass();
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    double prec_mass =
        ms_six_ptr_vec[0]->getMsHeaderPtr()->getFirstPrecMonoMassMinusWater();
    ProteoformPtr proteo_ptr = prsm_ptrs[i]->getProteoformPtr();
    int unexpect_shift_num = proteo_ptr->getAlterNum(AlterType::UNEXPECTED);
    ProteoformTypePtr type_ptr = proteo_ptr->getProteoformType();

    // the variable PTMs of the mass shifts, and the truncation N-terminal
    // modification, counted with its own type number
    int var_ptm_num = proteo_ptr->getAlterNum(AlterType::VARIABLE);
    int n_term_var_ptm_num = proteo_ptr->hasTruncNTermMod() ? 1 : 0;
    double search_mass = prec_mass;
    if (unexpect_shift_num == 0 && var_ptm_num == 0 &&
        n_term_var_ptm_num == 0) {
      // in ZERO PTM searching, +/-1 Da was allowed.
      // We need to adjust the prec mass for candidate number computation
      // if there was 1 Da difference between the precursor mass of the PRM
      // spectra and the adjusted prec mass (process passes the PrSMs of an
      // offset with the PRM spectra of the adjusted precursor, so that the
      // offset is not applied twice).
      if (std::abs(spec_prec_mass - prsm_ptrs[i]->getAdjustedPrecMass()) >
          tolerance) {
        if (spec_prec_mass < prsm_ptrs[i]->getAdjustedPrecMass()) {
          search_mass = prec_mass + mass_constant::getIsotopeMass();
        } else {
          search_mass = prec_mass - mass_constant::getIsotopeMass();
        }
      }
    }
    // when variable PTMs are allowed, use residue seq mass as precursor mass,
    // without the truncation N-terminal modification (the residue sequence
    // has it as its N-terminal modification)
    if (unexpect_shift_num == 0 && var_ptm_num + n_term_var_ptm_num > 0) {
      search_mass = proteo_ptr->getResSeqPtr()->getResMassSum() -
                    proteo_ptr->getTruncNTermModPtr()->getShift();
    }
    LOG_DEBUG("prec_ mass " << prec_mass << " search mass " << search_mass);
    double cand_num = test_num_ptr_->compCandNum(type_ptr, unexpect_shift_num,
                                                 search_mass, tolerance);
    if (cand_num == 0.0) {
      LOG_ERROR("Zero candidate number!");
      // cand_num = ExpectedValue::getMaxDouble();
      cand_num = 1.0;
    }

    // multiple a factor for variable PTMs: the truncation N-terminal
    // modification is counted with its own type number
    double var_ptm_factor =
        pow(mng_ptr_->var_ptm_type_num_ + 1, var_ptm_num) *
        pow(mng_ptr_->n_term_var_ptm_type_num_ + 1, n_term_var_ptm_num);
    LOG_DEBUG("candidate number "
              << cand_num << " var_ptm_type_num " << mng_ptr_->var_ptm_type_num_
              << " var ptm num " << var_ptm_num << " factor "
              << var_ptm_factor);
    cand_num = cand_num * var_ptm_factor;

    double adjust_factor = 1.0;
    if (type_ptr == ProteoformType::COMPLETE ||
        type_ptr == ProteoformType::PREFIX) {
      ExpectedValuePtr ev_ptr = std::make_shared<ExpectedValue>(
          prot_probs[i], cand_num, adjust_factor);
      prsm_ptrs[i]->setExpectedValuePtr(ev_ptr);
    } else {
      ExpectedValuePtr ev_ptr = std::make_shared<ExpectedValue>(
          pep_probs[i], cand_num, adjust_factor);
      prsm_ptrs[i]->setExpectedValuePtr(ev_ptr);
    }
  }
}

void CompPValueArray::compSingleExpectedValue(const DeconvMsPtrVec& ms_ptr_vec,
                                              const PrsmPtr& prsm_ptr,
                                              double n_term_label_mass) {
  // the adjusted precursor mass of a PrSM is label-free (toppic -N); the PRM
  // spectra are built from the observed one and subtract the label
  double refine_prec_mass = prsm_ptr->getAdjustedPrecMass() + n_term_label_mass;
  DeconvMsPtrVec refine_ms_ptr_vec =
      deconv_ms_util::getRefineMsPtrVec(ms_ptr_vec, refine_prec_mass);

  PrmMsPtrVec prm_ms_ptr_vec = prm_ms_factory::geneMsSixPtrVec(
      refine_ms_ptr_vec, mng_ptr_->prsm_para_ptr_->getSpParaPtr(),
      refine_prec_mass, n_term_label_mass);

  PrsmPtrVec prsm_ptrs;
  prsm_ptrs.push_back(prsm_ptr);
  compMultiExpectedValues(prm_ms_ptr_vec, prsm_ptrs, true,
                          useModResidues(prsm_ptr));
}

int CompPValueArray::getIsotopeOffset(const SpectrumSetPtr& spec_set_ptr,
                                      const PrsmPtr& prsm_ptr) {
  ProteoformPtr proteo_ptr = prsm_ptr->getProteoformPtr();
  if (proteo_ptr->getAlterNum(AlterType::UNEXPECTED) > 0 ||
      proteo_ptr->getAlterNum(AlterType::VARIABLE) > 0 ||
      proteo_ptr->hasTruncNTermMod()) {
    return 0;
  }
  // the rule of compMultiExpectedValues's candidate number: a PrSM without
  // mass shifts whose adjusted precursor mass differs from the spectrum's
  // (label-free) precursor mass by more than the precursor error tolerance
  // was matched at one isotope above or below it
  MsHeaderPtr header_ptr = spec_set_ptr->getMsSixPtrVec()[0]->getMsHeaderPtr();
  double prec_mass = header_ptr->getFirstPrecMonoMass();
  double tolerance = prec_tole_ptr_->compPrecErrorTole(prec_mass);
  double adjusted_mass = prsm_ptr->getAdjustedPrecMass();
  if (std::abs(prec_mass - adjusted_mass) <= tolerance) {
    return 0;
  }
  return prec_mass < adjusted_mass ? 1 : -1;
}

void CompPValueArray::process(const SpectrumSetPtr& spec_set_ptr,
                              PrsmPtrVec& prsm_ptrs, bool is_separate) {
  if (is_separate) {
    double n_term_label_mass = spec_set_ptr->getNTermLabelMass();
    for (unsigned i = 0; i < prsm_ptrs.size(); i++) {
      compSingleExpectedValue(spec_set_ptr->getDeconvMsPtrVec(), prsm_ptrs[i],
                              n_term_label_mass);
    }
    return;
  }
  // The zero-shift search scores a spectrum at its precursor mass P and at
  // P +- 1 isotope (SpPara::getZeroShiftSearchPrecErrorVec). A PrSM found
  // at P +- 1 gets its E-value on the PRM spectra of the spectrum set the
  // search built for that mass, with that mass as the generating function's
  // precursor mass, as if the precursor were reported there; the other
  // PrSMs share the generating function of the spectrum set's PRM spectra
  // (toppic_test experiment 117).
  std::map<int, PrsmPtrVec> offset_prsms;
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    int offset = getIsotopeOffset(spec_set_ptr, prsm_ptrs[i]);
    offset_prsms[offset].push_back(prsm_ptrs[i]);
  }
  for (auto& offset_prsm : offset_prsms) {
    int offset = offset_prsm.first;
    if (offset == 0) {
      // with the variable PTMs of var_ptm_mod_ptr_vec_, the variable PTM
      // PrSMs (all of offset 0) get the generating function of the
      // residues with their modified forms, the others main's (toppic_test
      // experiment 129)
      PrsmPtrVec var_prsm_ptrs;
      PrsmPtrVec other_prsm_ptrs;
      for (size_t i = 0; i < offset_prsm.second.size(); i++) {
        if (useModResidues(offset_prsm.second[i])) {
          var_prsm_ptrs.push_back(offset_prsm.second[i]);
        } else {
          other_prsm_ptrs.push_back(offset_prsm.second[i]);
        }
      }
      if (!other_prsm_ptrs.empty()) {
        compMultiExpectedValues(spec_set_ptr->getMsSixPtrVec(),
                                other_prsm_ptrs, false);
      }
      if (!var_prsm_ptrs.empty()) {
        compMultiExpectedValues(spec_set_ptr->getMsSixPtrVec(), var_prsm_ptrs,
                                false, true);
      }
      continue;
    }
    // the PRM spectra of spectrum_set_factory::geneSpectrumSetPtr for the
    // observed precursor mass + offset isotopes, as the search built them
    double prec_mono_mass = spec_set_ptr->getPrecMonoMass() +
                            offset * mass_constant::getIsotopeMass();
    PrmMsPtrVec ms_six_ptr_vec = prm_ms_factory::geneMsSixPtrVec(
        spec_set_ptr->getDeconvMsPtrVec(),
        mng_ptr_->prsm_para_ptr_->getSpParaPtr(), prec_mono_mass,
        spec_set_ptr->getNTermLabelMass());
    compMultiExpectedValues(ms_six_ptr_vec, offset_prsm.second, false);
  }
}

}  // namespace toppic
