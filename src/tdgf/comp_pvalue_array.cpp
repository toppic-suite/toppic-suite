#include "base/logger.hpp"
#include "tdgf/comp_pvalue_array.hpp"

namespace prot {

CompPValueArray::CompPValueArray(const ProteoformPtrVec &raw_proteo_ptrs, 
                                 const ProteoformPtrVec &mod_proteo_ptrs,
                                 const ResFreqPtrVec &residue_ptrs,
                                 TdgfMngPtr mng_ptr) {
  mng_ptr_ = mng_ptr;
  pep_n_term_residue_ptrs_ = residue_ptrs;
  prot_n_term_residue_ptrs_ = compNTermResidueFreq(mod_proteo_ptrs);
  LOG_DEBUG("n term residue_ptrs initialized")
  comp_prob_ptr_ = CompProbValuePtr(
      new CompProbValue(mng_ptr_->convert_ratio_,
                        residue_ptrs, mng_ptr_->unexpected_shift_num_ + 1, 
                        mng_ptr_->max_table_height_, 
                        mng_ptr_->max_prec_mass_));
  LOG_DEBUG("comp prob value initialized")
                        
  test_num_ptr_ = CountTestNumPtr(new CountTestNum(raw_proteo_ptrs, mod_proteo_ptrs,
                                                   residue_ptrs, mng_ptr->convert_ratio_,
                                                   mng_ptr->max_prec_mass_,
                                                   mng_ptr->max_ptm_mass_));
  LOG_DEBUG("test number initialized")
}

/* set alignment */
ExtremeValuePtrVec CompPValueArray::compExtremeValues(PrmMsPtr ms_six_ptr, 
                                                      const PrsmPtrVec &prsm_ptrs, 
                                                      bool strict) {
  PrmPeakPtrVec prm_peak_ptrs = ms_six_ptr->getPeakPtrVec();
  std::vector<double> prot_probs; 
  compProbArray(comp_prob_ptr_, prot_n_term_residue_ptrs_, 
                prm_peak_ptrs, prsm_ptrs, strict, prot_probs);
  std::vector<double> pep_probs;
  compProbArray(comp_prob_ptr_, pep_n_term_residue_ptrs_, 
                prm_peak_ptrs, prsm_ptrs, strict, pep_probs);
  //LOG_DEBUG("probability computation complete");
  double prec_mass = ms_six_ptr->getHeaderPtr()->getPrecMonoMassMinusWater();
  double tolerance = ms_six_ptr->getHeaderPtr()->getErrorTolerance();
  ExtremeValuePtrVec ev_prob_ptrs; 
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    //LOG_DEBUG("prsm " << i << " prsm size " << prsm_ptrs.size());
    int unexpect_shift_num = prsm_ptrs[i]->getProteoformPtr()->getUnexpectedChangeNum();
    SemiAlignTypePtr type_ptr = prsm_ptrs[i]->getProteoformPtr()->getSemiAlignType();
    double cand_num = test_num_ptr_->compCandNum(type_ptr, unexpect_shift_num, 
                                                 prec_mass, tolerance);
    //LOG_DEBUG("candidate number " << cand_num);
    if (cand_num == 0.0) {
      LOG_WARN("Zero candidate number");
      cand_num = std::numeric_limits<double>::infinity();
    }
    if (type_ptr == SemiAlignTypeFactory::getCompletePtr() 
        || type_ptr == SemiAlignTypeFactory::getPrefixPtr()) {
      ev_prob_ptrs.push_back(ExtremeValuePtr(new ExtremeValue(prot_probs[i], cand_num, 1)));
    } else {
      ev_prob_ptrs.push_back(ExtremeValuePtr(new ExtremeValue(pep_probs[i], cand_num,1)));
    }
    //LOG_DEBUG("assignment complete");
  }
  return ev_prob_ptrs;
}

ExtremeValuePtr CompPValueArray::compExtremeValue(PrmMsPtr prm_ms_ptr, 
                                                  PrsmPtr prsm_ptr) {
  PrsmPtrVec prsm_ptrs;
  prsm_ptrs.push_back(prsm_ptr);
  ExtremeValuePtrVec extreme_values = compExtremeValues(prm_ms_ptr, prsm_ptrs, true);
  return extreme_values[0];
}

void CompPValueArray::setPValue(DeconvMsPtr ms_ptr, PrsmPtr prsm_ptr) {
  double refine_prec_mass = prsm_ptr->getAdjustedPrecMass();
  DeconvMsPtr refine_ms_ptr = getRefineMs(ms_ptr, prsm_ptr->getCalibration(),
                                          refine_prec_mass);

  LOG_DEBUG("recalibration " << prsm_ptr->getCalibration()
               << " original precursor "
               << ms_ptr->getHeaderPtr()->getPrecMonoMass()
               << " precursor " << refine_prec_mass);
  // Delta = 0 is important 
  PrmMsPtr prm_ms_ptr = createMsSixPtr(refine_ms_ptr, 0, mng_ptr_->prsm_para_ptr_->getSpParaPtr());
  ExtremeValuePtr prob_ptr = compExtremeValue(prm_ms_ptr, prsm_ptr);
  prsm_ptr->setProbPtr(prob_ptr);
}

void CompPValueArray::setPValueArray(PrmMsPtr prm_ms_ptr, PrsmPtrVec &prsm_ptrs) {
  ExtremeValuePtrVec extreme_values = compExtremeValues(prm_ms_ptr, prsm_ptrs, false);
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    prsm_ptrs[i]->setProbPtr(extreme_values[i]);
  }
  //LOG_DEBUG("Set value complete");
}

}
