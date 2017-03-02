// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include "base/logger.hpp"
#include "spec/deconv_ms_factory.hpp"
#include "tdgf/comp_pvalue_array.hpp"

namespace prot {

CompPValueArray::CompPValueArray(CountTestNumPtr test_num_ptr,
                                 TdgfMngPtr mng_ptr) {
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
  //LOG_DEBUG("probability computation complete");
  double prec_mass = ms_six_ptr_vec[0]->getMsHeaderPtr()->getPrecMonoMassMinusWater();
  double tolerance = ms_six_ptr_vec[0]->getMsHeaderPtr()->getErrorTolerance(ppo);
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    //LOG_DEBUG("prsm " << i << " prsm size " << prsm_ptrs.size());
    int unexpect_shift_num = prsm_ptrs[i]->getProteoformPtr()->getChangeNum(ChangeType::UNEXPECTED);
    AlignTypePtr type_ptr = prsm_ptrs[i]->getProteoformPtr()->getAlignType();

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

    double cand_num = test_num_ptr_->compCandNum(type_ptr, index, prec_mass, tolerance);

    //LOG_DEBUG("Shift number " << unexpect_shift_num << " type " << type_ptr->getName() 
    //<< " one prob " << prot_probs[i] << " cand num " << cand_num);
    //LOG_DEBUG("candidate number " << cand_num);
    if (cand_num == 0.0) {
      LOG_WARN("Zero candidate number");
      cand_num = std::numeric_limits<double>::infinity();
    }

    if (type_ptr == AlignType::COMPLETE || type_ptr == AlignType::PREFIX) {
      ExtremeValuePtr ev_ptr = std::make_shared<ExtremeValue>(prot_probs[i], cand_num, 1);
      prsm_ptrs[i]->setExtremeValuePtr(ev_ptr);
    } else {
      ExtremeValuePtr ev_ptr = std::make_shared<ExtremeValue>(pep_probs[i], cand_num, 1);
      prsm_ptrs[i]->setExtremeValuePtr(ev_ptr);
    }
    //LOG_DEBUG("assignment complete");
  }
}

void CompPValueArray::compSingleExtremeValue(const DeconvMsPtrVec &ms_ptr_vec, 
                                             PrsmPtr prsm_ptr, double ppo) {
  double refine_prec_mass = prsm_ptr->getAdjustedPrecMass();
  DeconvMsPtrVec refine_ms_ptr_vec = DeconvMsFactory::getRefineMsPtrVec(ms_ptr_vec, refine_prec_mass);

  /*
     LOG_DEBUG("recalibration " << prsm_ptr->getCalibration()
     << " original precursor "
     << ms_ptr->getHeaderPtr()->getPrecMonoMass()
     << " precursor " << refine_prec_mass);
     */
  PrmMsPtrVec prm_ms_ptr_vec = PrmMsFactory::geneMsSixPtrVec(refine_ms_ptr_vec, 
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
      compSingleExtremeValue(spec_set_ptr->getDeconvMsPtrVec(), 
                             prsm_ptrs[i], ppo);
    }
  } else {
    compMultiExtremeValues(spec_set_ptr->getMsSixPtrVec(), prsm_ptrs, ppo, false);
  }
}

}
