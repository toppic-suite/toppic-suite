#include "base/logger.hpp"
#include "tdgf/comp_pvalue_array.hpp"

namespace prot {

CompPValueArray::CompPValueArray(ProteoformPtrVec &raw_forms, 
                                 ProteoformPtrVec &prot_mod_forms,
                                 ResFreqPtrVec &residues,
                                 TdgfMngPtr mng_ptr) {
  mng_ptr_ = mng_ptr;
  pep_n_term_residues_ = residues;
  prot_n_term_residues_ = compNTermResidueFreq(prot_mod_forms);
  LOG_DEBUG("n term residues initialized")
  comp_prob_ptr_ = CompProbValuePtr(
      new CompProbValue(mng_ptr_->double_to_int_constant_,
                        residues, mng_ptr_->unexpected_shift_num_ + 1, 
                        mng_ptr_->max_table_height_, 
                        mng_ptr_->max_sp_prec_mass_));
  LOG_DEBUG("comp prob value initialized")
                        
  test_num_ptr_ = CountTestNumPtr(new CountTestNum(raw_forms, prot_mod_forms,
                                                   residues, mng_ptr));
  LOG_DEBUG("test number initialized")
}

/* set alignment */
ExtremeValuePtrVec CompPValueArray::compExtremeValues(PrmMsPtr ms_six, 
                                                      PrSMPtrVec &prsms, 
                                                      bool strict) {
  PrmPeakPtrVec prm_peaks = ms_six->getPeakPtrVec();
  std::vector<double> prot_probs; 
  compProbArray(comp_prob_ptr_, prot_n_term_residues_, 
                prm_peaks, prsms, strict, prot_probs);
  std::vector<double> pep_probs;
  compProbArray(comp_prob_ptr_, pep_n_term_residues_, 
                prm_peaks, prsms, strict, pep_probs);
  //LOG_DEBUG("probability computation complete");
  double prec_mass = ms_six->getHeaderPtr()->getPrecMonoMassMinusWater();
  double tolerance = ms_six->getHeaderPtr()->getErrorTolerance();
  ExtremeValuePtrVec ev_probs; 
  for (unsigned int i = 0; i < prsms.size(); i++) {
    //LOG_DEBUG("prsm " << i << " prsm size " << prsms.size());
    int unexpect_shift_num = prsms[i]->getProteoformPtr()->getUnexpectedChangeNum();
    SemiAlignTypePtr t = prsms[i]->getProteoformPtr()->getSemiAlignType();
    double cand_num = test_num_ptr_->compCandNum(t, unexpect_shift_num, 
                                                 prec_mass, tolerance);
    //LOG_DEBUG("candidate number " << cand_num);
    if (cand_num == 0.0) {
      LOG_WARN("Zero candidate number");
      cand_num = std::numeric_limits<double>::infinity();
    }
    if (t == SemiAlignTypeFactory::getCompletePtr() 
        || t == SemiAlignTypeFactory::getPrefixPtr()) {
      ev_probs.push_back(ExtremeValuePtr(new ExtremeValue(prot_probs[i], cand_num, 1)));
      /*
      logger.debug("AlignType " + t + " shift "
                   + prsms[i].getAnnoProtein().getUnknownShiftNum()
                   + " thresh " + prsms[i].getMatchFragNum()
                   + " prot one prob " + protProbs[i] + " pep one prob "
                   + pepProbs[i] + " candiate num " + nCandidates
                   + " p value " + eVProbs[i].getPValue());
                   */
    } else {
      ev_probs.push_back(ExtremeValuePtr(new ExtremeValue(pep_probs[i], cand_num,1)));
      /*
      logger.debug("AlignType " + t + " shift "
                   + prsms[i].getAnnoProtein().getUnknownShiftNum()
                   + " thresh " + prsms[i].getMatchFragNum()
                   + " one prob " + pepProbs[i] + " candiate num "
                   + nCandidates + " p value " + eVProbs[i].getPValue());
                   */
    }
    //LOG_DEBUG("assignment complete");
  }
  return ev_probs;
}

ExtremeValuePtr CompPValueArray::compExtremeValue(PrmMsPtr prm_ms_ptr, 
                                                  PrSMPtr prsm_ptr) {
  PrSMPtrVec prsms;
  prsms.push_back(prsm_ptr);
  ExtremeValuePtrVec extreme_values = compExtremeValues(prm_ms_ptr, prsms, true);
  return extreme_values[0];
}

void CompPValueArray::setPValue(DeconvMsPtr ms_ptr, PrSMPtr prsm_ptr) {
  double refine_prec_mass = prsm_ptr->getAdjustedPrecMass();
  DeconvMsPtr refine_ms_ptr = getRefineMs(ms_ptr, prsm_ptr->getCalibration(),
                                          refine_prec_mass);

  LOG_DEBUG("recalibration " << prsm_ptr->getCalibration()
               << " original precursor "
               << ms_ptr->getHeaderPtr()->getPrecMonoMass()
               << " precursor " << refine_prec_mass);
  // Delta = 0 is important 
  PrmMsPtr prm_ms_ptr = getMsSix(refine_ms_ptr, 0, mng_ptr_->prsm_para_ptr_->getSpParaPtr());
  ExtremeValuePtr prob_ptr = compExtremeValue(prm_ms_ptr, prsm_ptr);
  prsm_ptr->setProbPtr(prob_ptr);
}

void CompPValueArray::setPValueArray(PrmMsPtr prm_ms_ptr, PrSMPtrVec prsms) {
  ExtremeValuePtrVec extreme_values = compExtremeValues(prm_ms_ptr, prsms, false);
  for (unsigned int i = 0; i < prsms.size(); i++) {
    prsms[i]->setProbPtr(extreme_values[i]);
  }
  //LOG_DEBUG("Set value complete");
}

}
