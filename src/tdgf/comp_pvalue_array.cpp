//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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
#include <limits>
#include <vector>

#include "util/logger.hpp"
#include "base/base_data.hpp"
#include "spec/deconv_ms_util.hpp"
#include "spec/prm_ms_factory.hpp"
#include "tdgf/comp_pvalue_array.hpp"

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
  LOG_DEBUG("comp prob value initialized")
}

/* set alignment */
void CompPValueArray::compMultiExtremeValues(const PrmMsPtrVec &ms_six_ptr_vec,
                                             PrsmPtrVec &prsm_ptrs,
                                             double ppo, bool strict) {
  PrmPeakPtrVec2D prm_ptr_2d;
  for (size_t i = 0; i < ms_six_ptr_vec.size(); i++) {
    PrmPeakPtrVec prm_peak_ptrs = ms_six_ptr_vec[i]->getPeakPtrVec();
    prm_ptr_2d.push_back(prm_peak_ptrs);
  }
  std::vector<double> prot_probs;
  CompProbValue::compProbArray(comp_prob_ptr_, prot_n_term_residue_ptrs_,
                               prm_ptr_2d, prsm_ptrs, strict, prot_probs);
  std::vector<double> pep_probs;
  CompProbValue::compProbArray(comp_prob_ptr_, pep_n_term_residue_ptrs_,
                               prm_ptr_2d, prsm_ptrs, strict, pep_probs);
  // LOG_DEBUG("probability computation complete");
  double tolerance = ms_six_ptr_vec[0]->getMsHeaderPtr()->getErrorTolerance(ppo);
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    double prec_mass = ms_six_ptr_vec[0]->getMsHeaderPtr()->getPrecMonoMassMinusWater();
    // LOG_DEBUG("prsm " << i << " prsm size " << prsm_ptrs.size());
    int unexpect_shift_num = prsm_ptrs[i]->getProteoformPtr()->getMassShiftNum(MassShiftType::UNEXPECTED);
    ProteoformTypePtr type_ptr = prsm_ptrs[i]->getProteoformPtr()->getProteoformType();

    if (unexpect_shift_num == 0) {
      // in ZERO PTM searching, +/-1 Da was allowed.
      // We need to adjust the prec mass for candidate number computation
      // if there was 1 Da difference between original prec mass and adjusted
      // prec mass.
      if (std::abs(prsm_ptrs[i]->getOriPrecMass() - prsm_ptrs[i]->getAdjustedPrecMass()) > tolerance) {
        if (prsm_ptrs[i]->getOriPrecMass() < prsm_ptrs[i]->getAdjustedPrecMass()) {
          prec_mass += 1;
        } else {
          prec_mass -= 1;
        }
      }
    }

    int index;
    if (mng_ptr_->variable_ptm_) {
      // need to be improved
      if (unexpect_shift_num == 0) {
        index = 1;
      } else {
        index = unexpect_shift_num;
      }
    } else {
      index = unexpect_shift_num;
    }

    LOG_DEBUG("prec_mass " << prec_mass);
    double cand_num = test_num_ptr_->compCandNum(type_ptr, index, prec_mass, tolerance);

    // LOG_DEBUG("Shift number " << unexpect_shift_num << " type " << type_ptr->getName()
    // << " one prob " << prot_probs[i] << " cand num " << cand_num);
    // LOG_DEBUG("candidate number " << cand_num);
    if (cand_num == 0.0) {
      LOG_WARN("Zero candidate number");
      cand_num = ExtremeValue::getMaxDouble();
    }

    if (type_ptr == ProteoformType::COMPLETE || type_ptr == ProteoformType::PREFIX) {
      ExtremeValuePtr ev_ptr = std::make_shared<ExtremeValue>(prot_probs[i], cand_num, 1);
      prsm_ptrs[i]->setExtremeValuePtr(ev_ptr);
    } else {
      ExtremeValuePtr ev_ptr = std::make_shared<ExtremeValue>(pep_probs[i], cand_num, 1);
      prsm_ptrs[i]->setExtremeValuePtr(ev_ptr);
    }
    // LOG_DEBUG("assignment complete");
  }
}

void CompPValueArray::compSingleExtremeValue(const DeconvMsPtrVec &ms_ptr_vec,
                                             PrsmPtr prsm_ptr, double ppo) {
  double refine_prec_mass = prsm_ptr->getAdjustedPrecMass();
  DeconvMsPtrVec refine_ms_ptr_vec = deconv_ms_util::getRefineMsPtrVec(ms_ptr_vec, refine_prec_mass);

  /*
     LOG_DEBUG("recalibration " << prsm_ptr->getCalibration()
     << " original precursor "
     << ms_ptr->getHeaderPtr()->getPrecMonoMass()
     << " precursor " << refine_prec_mass);
     */
  PrmMsPtrVec prm_ms_ptr_vec = prm_ms_factory::geneMsSixPtrVec(refine_ms_ptr_vec,
                                                               mng_ptr_->prsm_para_ptr_->getSpParaPtr(),
                                                               refine_prec_mass);

  PrsmPtrVec prsm_ptrs;
  prsm_ptrs.push_back(prsm_ptr);
  compMultiExtremeValues(prm_ms_ptr_vec, prsm_ptrs, ppo, true);
}

void CompPValueArray::process(SpectrumSetPtr spec_set_ptr, PrsmPtrVec &prsm_ptrs,
                              double ppo, bool is_separate) {
  if (is_separate) {
    for (unsigned i = 0; i < prsm_ptrs.size(); i++) {
      compSingleExtremeValue(spec_set_ptr->getDeconvMsPtrVec(), prsm_ptrs[i], ppo);
    }
  } else {
    compMultiExtremeValues(spec_set_ptr->getMsSixPtrVec(), prsm_ptrs, ppo, false);
  }
}

}  // namespace toppic
