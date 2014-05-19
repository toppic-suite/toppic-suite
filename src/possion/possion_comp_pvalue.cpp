#include "base/logger.hpp"
#include "tdgf/comp_prob_value.hpp"
#include "possion/possion_comp_pvalue.hpp"

namespace prot {

PossionCompPValue::PossionCompPValue(ProteoformPtrVec &raw_forms, 
                                     ProteoformPtrVec &prot_mod_forms,
                                     PossionMngPtr mng_ptr) {
  mng_ptr_ = mng_ptr;
  ResFreqPtrVec residue_freqs 
      = compResidueFreq(mng_ptr_->prsm_para_ptr_->getFixModResiduePtrVec(), raw_forms); 
  test_num_ptr_ = CountTestNumPtr(new CountTestNum(raw_forms, prot_mod_forms,
                                                   residue_freqs,
                                                   mng_ptr->convert_ratio_,
                                                   mng_ptr->max_prec_mass_,
                                                   mng_ptr->max_ptm_mass_));

  residue_avg_len_ = computeAvgLength(residue_freqs, 1.0);
  LOG_DEBUG("test number initialized")
}

double PossionCompPValue::compRandMatchProb(double prec_mass, bool is_strict) {
  double avg_mass = prec_mass / 2.0;
  PeakTolerancePtr tole_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getPeakTolerancePtr();
  // error tolerance for N-term + error tolerance for C -term
  double coverage;
  if (is_strict) {
    coverage = tole_ptr->compStrictErrorTole(avg_mass) * 2 + tole_ptr->compStrictErrorTole(avg_mass) * 2;
  }
  else {
    coverage = tole_ptr->compStrictErrorTole(avg_mass) * 2 + tole_ptr->compRelaxErrorTole(prec_mass, avg_mass) * 2;
  }
  coverage = coverage * mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getExtendOffsets().size();
  double prob = coverage / residue_avg_len_;
  std::cout << " Random match prob " << prob << std::endl;
  return prob;
}

double PossionCompPValue::compConditionProb(double rand_match_prob, PrmMsPtr ms_six, PrsmPtr prsm) {
  double f = ms_six->size();
  double x = rand_match_prob;
  double n = prsm->getMatchFragNum();
  n = n - prsm->getProteoformPtr()->getUnexpectedChangeNum() * mng_ptr_->shift_penalty;
  if (f <= 0.0 || x <= 0.0 || n <= 0.0) {
    return 1;
  }
  double xf = x*f;
  double sum = 0;
  double xf_power = 1;
  double factorial = 1;
  for (int i = 0; i < n; i++) {
    sum = sum + std::exp(-xf) * xf_power / factorial;
    xf_power = xf_power * xf;
    factorial = factorial * (i+1);
  }
  double prob = 1 - sum;
  std::cout << " f " << f << " x " << x  << " n " << n << " prob " << prob  << std::endl;
  return prob;
}

ExtremeValuePtrVec PossionCompPValue::compExtremeValues(PrmMsPtr ms_six, 
                                                        PrsmPtrVec &prsms,
                                                        bool is_strict) {
  double prec_mass = ms_six->getHeaderPtr()->getPrecMonoMassMinusWater();
  double tolerance = ms_six->getHeaderPtr()->getErrorTolerance();
  double rand_match_prob = compRandMatchProb(prec_mass, is_strict); 

  ExtremeValuePtrVec ev_probs; 
  for (unsigned int i = 0; i < prsms.size(); i++) {
    double cond_prob = compConditionProb(rand_match_prob, ms_six, prsms[i]);
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
      ev_probs.push_back(ExtremeValuePtr(new ExtremeValue(cond_prob, cand_num, 1)));
    } else {
      ev_probs.push_back(ExtremeValuePtr(new ExtremeValue(cond_prob, cand_num,1)));
    }
    //LOG_DEBUG("assignment complete");
  }
  return ev_probs;
}

void PossionCompPValue::setPValueArray(PrmMsPtr prm_ms_ptr, PrsmPtrVec &prsms) {
  ExtremeValuePtrVec extreme_values = compExtremeValues(prm_ms_ptr, prsms, false);
  for (unsigned int i = 0; i < prsms.size(); i++) {
    prsms[i]->setProbPtr(extreme_values[i]);
  }
  //LOG_DEBUG("Set value complete");
}

}
