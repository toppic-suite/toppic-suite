//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#include <cmath>

#include "common/util/logger.hpp"
#include "common/base/mass_constant.hpp"
#include "ms/spec/deconv_ms_util.hpp"
#include "ms/factory/prm_ms_factory.hpp"
#include "stat/tdgf/comp_prob_value_array.hpp"
#include "stat/tdgf/comp_pvalue_array.hpp"

namespace toppic {

CompPValueArray::CompPValueArray(CountTestNumPtr test_num_ptr, TdgfMngPtr mng_ptr) {
  mng_ptr_ = mng_ptr;
  test_num_ptr_ = test_num_ptr;
  residue_ptrs_ = test_num_ptr->getResFreqPtrVec();
  pep_n_term_residue_ptrs_ = residue_ptrs_;
  prot_n_term_residue_ptrs_ = test_num_ptr->getNTermResFreqPtrVec();
  comp_prob_ptr_ = std::make_shared<CompProbValue>(mng_ptr_->convert_ratio_,
                                                   residue_ptrs_,
                                                   test_num_ptr->getResidueAvgLen(),
                                                   mng_ptr_->unexpected_shift_num_ + 1,
                                                   mng_ptr_->max_table_height_,
                                                   mng_ptr_->max_prec_mass_);
}

// set alignment 
void CompPValueArray::compMultiExpectedValues(const PrmMsPtrVec &ms_six_ptr_vec,
                                              PrsmPtrVec &prsm_ptrs,
                                              double ppo, bool strict) {
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
  }
  else {
    prob_prec_mass = ms_six_ptr_vec[0]->getMsHeaderPtr()->getPrecMonoMass();
  }
  PeakTolerancePtr tole_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr();
  std::vector<double> prot_probs;
  comp_prob_value_array::compProbArray(comp_prob_ptr_, prot_n_term_residue_ptrs_,
                                      prm_ptr_2d, prsm_ptrs, strict, 
                                      prob_prec_mass, tole_ptr, prot_probs);
  std::vector<double> pep_probs;
  comp_prob_value_array::compProbArray(comp_prob_ptr_, pep_n_term_residue_ptrs_,
                               prm_ptr_2d, prsm_ptrs, strict, prob_prec_mass, tole_ptr, pep_probs);
  double tolerance = ms_six_ptr_vec[0]->getMsHeaderPtr()->getErrorTolerance(ppo);
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    double prec_mass = ms_six_ptr_vec[0]->getMsHeaderPtr()->getPrecMonoMassMinusWater();
    int unexpect_shift_num = prsm_ptrs[i]->getProteoformPtr()->getAlterNum(AlterType::UNEXPECTED);
    ProteoformTypePtr type_ptr = prsm_ptrs[i]->getProteoformPtr()->getProteoformType();

    if (unexpect_shift_num == 0) {
      // in ZERO PTM searching, +/-1 Da was allowed.
      // We need to adjust the prec mass for candidate number computation
      // if there was 1 Da difference between original prec mass and adjusted
      // prec mass.
      if (std::abs(prsm_ptrs[i]->getOriPrecMass() - prsm_ptrs[i]->getAdjustedPrecMass()) > tolerance) {
        if (prsm_ptrs[i]->getOriPrecMass() < prsm_ptrs[i]->getAdjustedPrecMass()) {
          prec_mass += mass_constant::getIsotopeMass();
        } 
        else {
          prec_mass -= mass_constant::getIsotopeMass();
        }
      }
    }

    int index;
    if (mng_ptr_->variable_ptm_) {
      // need to be improved
      if (unexpect_shift_num == 0) {
        index = 1;
      } 
      else {
        index = unexpect_shift_num;
      }
    } 
    else {
      index = unexpect_shift_num;
    }

    double cand_num = test_num_ptr_->compCandNum(type_ptr, index, prec_mass, tolerance);

    if (cand_num == 0.0) {
      LOG_WARN("Zero candidate number");
      cand_num = ExpectedValue::getMaxDouble();
    }

    if (type_ptr == ProteoformType::COMPLETE || type_ptr == ProteoformType::PREFIX) {
      ExpectedValuePtr ev_ptr = std::make_shared<ExpectedValue>(prot_probs[i], cand_num, 1);
      prsm_ptrs[i]->setExpectedValuePtr(ev_ptr);
    } 
    else {
      ExpectedValuePtr ev_ptr = std::make_shared<ExpectedValue>(pep_probs[i], cand_num, 1);
      prsm_ptrs[i]->setExpectedValuePtr(ev_ptr);
    }
  }
}

void CompPValueArray::compSingleExpectedValue(const DeconvMsPtrVec &ms_ptr_vec,
                                             PrsmPtr prsm_ptr, double ppo) {
  double refine_prec_mass = prsm_ptr->getAdjustedPrecMass();
  DeconvMsPtrVec refine_ms_ptr_vec = deconv_ms_util::getRefineMsPtrVec(ms_ptr_vec, refine_prec_mass);

  PrmMsPtrVec prm_ms_ptr_vec = prm_ms_factory::geneMsSixPtrVec(refine_ms_ptr_vec,
                                                               mng_ptr_->prsm_para_ptr_->getSpParaPtr(),
                                                               refine_prec_mass);

  PrsmPtrVec prsm_ptrs;
  prsm_ptrs.push_back(prsm_ptr);
  compMultiExpectedValues(prm_ms_ptr_vec, prsm_ptrs, ppo, true);
}

void CompPValueArray::process(SpectrumSetPtr spec_set_ptr, PrsmPtrVec &prsm_ptrs,
                              double ppo, bool is_separate) {
  if (is_separate) {
    for (unsigned i = 0; i < prsm_ptrs.size(); i++) {
      compSingleExpectedValue(spec_set_ptr->getDeconvMsPtrVec(), prsm_ptrs[i], ppo);
    }
  } 
  else {
    compMultiExpectedValues(spec_set_ptr->getMsSixPtrVec(), prsm_ptrs, ppo, false);
  }
}

}  // namespace toppic