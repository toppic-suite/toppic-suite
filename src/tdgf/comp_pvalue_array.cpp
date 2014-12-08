#include "base/logger.hpp"
#include "tdgf/comp_pvalue_array.hpp"

namespace prot {

CompPValueArray::CompPValueArray(CountTestNumPtr test_num_ptr,
                                 TdgfMngPtr mng_ptr) {
  mng_ptr_ = mng_ptr;
  test_num_ptr_ = test_num_ptr;
  residue_ptrs_ = test_num_ptr->getResFreqPtrVec();
  pep_n_term_residue_ptrs_ = residue_ptrs_;
  prot_n_term_residue_ptrs_ = test_num_ptr->getNTermResFreqPtrVec();
  comp_prob_ptr_ = CompProbValuePtr(
      new CompProbValue(mng_ptr_->convert_ratio_,
                        residue_ptrs_, 
                        test_num_ptr->getResidueAvgLen(),
                        mng_ptr_->unexpected_shift_num_ + 1, 
                        mng_ptr_->max_table_height_, 
                        mng_ptr_->max_prec_mass_));
  LOG_DEBUG("comp prob value initialized")
}

/* set alignment */
void CompPValueArray::compMultiExtremeValues(const PrmMsPtrVec &ms_six_ptr_vec, 
                                             PrsmPtrVec &prsm_ptrs, 
                                             bool strict) {
  PrmPeakPtrVec2D prm_ptr_2d;
  for (size_t i = 0; i < ms_six_ptr_vec.size(); i++) {
    PrmPeakPtrVec prm_peak_ptrs = ms_six_ptr_vec[i]->getPeakPtrVec();
    prm_ptr_2d.push_back(prm_peak_ptrs);
  }
  std::vector<double> prot_probs;
  compProbArray(comp_prob_ptr_, prot_n_term_residue_ptrs_, 
                prm_ptr_2d, prsm_ptrs, strict, prot_probs);
  std::vector<double> pep_probs;
  compProbArray(comp_prob_ptr_, pep_n_term_residue_ptrs_, 
                prm_ptr_2d, prsm_ptrs, strict, pep_probs);
  //LOG_DEBUG("probability computation complete");
  double prec_mass = ms_six_ptr_vec[0]->getHeaderPtr()->getPrecMonoMassMinusWater();
  double tolerance = ms_six_ptr_vec[0]->getHeaderPtr()->getErrorTolerance();
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    //LOG_DEBUG("prsm " << i << " prsm size " << prsm_ptrs.size());
    int unexpect_shift_num = prsm_ptrs[i]->getProteoformPtr()->getUnexpectedChangeNum();
    SemiAlignTypePtr type_ptr = prsm_ptrs[i]->getProteoformPtr()->getSemiAlignType();
    double cand_num = test_num_ptr_->compCandNum(type_ptr, unexpect_shift_num, 
                                                 prec_mass, tolerance);
    //LOG_DEBUG("Shift number " << unexpect_shift_num << " type " << type_ptr->getName() 
    //<< " one prob " << prot_probs[i] << " cand num " << cand_num);
    //LOG_DEBUG("candidate number " << cand_num);
    if (cand_num == 0.0) {
      LOG_WARN("Zero candidate number");
      cand_num = std::numeric_limits<double>::infinity();
    }
    if (type_ptr == SemiAlignTypeFactory::getCompletePtr() 
        || type_ptr == SemiAlignTypeFactory::getPrefixPtr()) {
      ExtremeValuePtr prob_ptr(new ExtremeValue(prot_probs[i], cand_num, 1));
      prsm_ptrs[i]->setProbPtr(prob_ptr);
    } else {
      ExtremeValuePtr prob_ptr(new ExtremeValue(pep_probs[i], cand_num, 1));
      prsm_ptrs[i]->setProbPtr(prob_ptr);
    }
    //LOG_DEBUG("assignment complete");
  }
}

void CompPValueArray::compSingleExtremeValue(const DeconvMsPtrVec &ms_ptr_vec, 
                                             PrsmPtr prsm_ptr) {
  double refine_prec_mass = prsm_ptr->getAdjustedPrecMass();
  DeconvMsPtrVec refine_ms_ptr_vec = getRefineMsPtrVec(ms_ptr_vec, refine_prec_mass);

  /*
  LOG_DEBUG("recalibration " << prsm_ptr->getCalibration()
               << " original precursor "
               << ms_ptr->getHeaderPtr()->getPrecMonoMass()
               << " precursor " << refine_prec_mass);
               */
  PrmMsPtrVec prm_ms_ptr_vec = createMsSixPtrVec(refine_ms_ptr_vec, 
                                                 mng_ptr_->prsm_para_ptr_->getSpParaPtr(),
                                                 refine_prec_mass);

  PrsmPtrVec prsm_ptrs;
  prsm_ptrs.push_back(prsm_ptr);
  compMultiExtremeValues(prm_ms_ptr_vec, prsm_ptrs, true);
}


void CompPValueArray::process(SpectrumSetPtr spec_set_ptr, 
                              bool is_separate, PrsmPtrVec &prsm_ptrs) {

  if (is_separate) {
    for (unsigned i = 0; i < prsm_ptrs.size(); i++) {
      compSingleExtremeValue(spec_set_ptr->getDeconvMsPtrVec(), 
                             prsm_ptrs[i]);
    }
  } 
  else {
    compMultiExtremeValues(spec_set_ptr->getMsSixPtrVec(), prsm_ptrs, false);
  }
}

}
