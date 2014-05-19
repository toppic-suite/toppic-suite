#include "base/logger.hpp"
#include "tdgf/comp_pvalue_poission.hpp"

namespace prot {

CompPValuePoission::CompPValuePoission(ProteoformPtrVec &raw_forms, 
                                       ProteoformPtrVec &prot_mod_forms,
                                       TdgfMngPtr mng_ptr) {
  mng_ptr_ = mng_ptr;
  test_num_ptr_ = CountTestNumPtr(new CountTestNum(raw_forms, prot_mod_forms,
                                                   residues, mng_ptr));
  LOG_DEBUG("test number initialized")
}

/* set alignment */
ExtremeValuePtrVec CompPValuePossion::compExtremeValues(PrmMsPtr ms_six, 
                                                        PrsmPtrVec &prsms) {
  PrmPeakPtrVec prm_peaks = ms_six->getPeakPtrVec();
  std::vector<double> prot_probs; 
  compProbArray(comp_prob_ptr_, prot_n_term_residues_, 
                prm_peaks, prsms, strict, prot_probs);
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
    } else {
      ev_probs.push_back(ExtremeValuePtr(new ExtremeValue(pep_probs[i], cand_num,1)));
    }
    //LOG_DEBUG("assignment complete");
  }
  return ev_probs;
}

ExtremeValuePtr CompPValueArray::compExtremeValue(PrmMsPtr prm_ms_ptr, 
                                                  PrsmPtr prsm_ptr) {
  PrsmPtrVec prsms;
  prsms.push_back(prsm_ptr);
  ExtremeValuePtrVec extreme_values = compExtremeValues(prm_ms_ptr, prsms, true);
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
  PrmMsPtr prm_ms_ptr = getMsSix(refine_ms_ptr, 0, mng_ptr_->prsm_para_ptr_->getSpParaPtr());
  ExtremeValuePtr prob_ptr = compExtremeValue(prm_ms_ptr, prsm_ptr);
  prsm_ptr->setProbPtr(prob_ptr);
}

void CompPValueArray::setPValueArray(PrmMsPtr prm_ms_ptr, PrsmPtrVec prsms) {
  ExtremeValuePtrVec extreme_values = compExtremeValues(prm_ms_ptr, prsms, false);
  for (unsigned int i = 0; i < prsms.size(); i++) {
    prsms[i]->setProbPtr(extreme_values[i]);
  }
  //LOG_DEBUG("Set value complete");
}

}
