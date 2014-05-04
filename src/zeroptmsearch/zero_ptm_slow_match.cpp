#include "base/logger.hpp"
#include "base/proteoform.hpp"
#include "base/activation.hpp"
#include "base/algorithm.hpp"
#include "spec/theo_peak.hpp"
#include "zeroptmsearch/zero_ptm_slow_match.hpp"

namespace prot {

ZeroPtmSlowMatch::ZeroPtmSlowMatch(DeconvMsPtr deconv_ms_ptr, 
                                   ZpFastMatchPtr fast_match_ptr,
                                   ZeroPtmMngPtr mng_ptr) {
  mng_ptr_ = mng_ptr;
  deconv_ms_ptr_ = deconv_ms_ptr;
  fast_match_ptr_ = fast_match_ptr;

  proteoform_ptr_ = getSubProteoform(fast_match_ptr->getProteoformPtr(), 
                                     fast_match_ptr->getBegin(), 
                                     fast_match_ptr->getEnd());

  SpParaPtr sp_para_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr();
  refine_prec_mass_ = proteoform_ptr_->getResSeqPtr()->getSeqMass();
  double delta = refine_prec_mass_ - deconv_ms_ptr->getHeaderPtr()->getPrecMonoMass();
  refine_ms_ptr_ = getMsThree(deconv_ms_ptr_, delta, sp_para_ptr);

  ActivationPtr activation_ptr = deconv_ms_ptr_->getHeaderPtr()->getActivationPtr();
  double min_mass = sp_para_ptr->getMinMass();
  TheoPeakPtrVec theo_peaks = getProteoformTheoPeak(proteoform_ptr_, 
                                                    activation_ptr, min_mass);

  compScore(refine_ms_ptr_,theo_peaks, sp_para_ptr->getPeakTolerancePtr()->getPpo());
}

// compute the average ppo
double compAvg(std::vector<double> ppos, double recal_ppo) {
  int cnt = 0;
  double sum = 0;
  for (unsigned int i = 0; i < ppos.size(); i++) {
    if (std::abs(ppos[i]) <= recal_ppo) {
      cnt++;
      sum += ppos[i];
    }
  }
  if (cnt == 0) {
    return 0;
  } else {
    return sum / cnt;
  }
}

/**
 * compute the validation of the candidates based on the different of
 * precursor mass
 */
bool ZeroPtmSlowMatch::isValid (double recal, double ppo) {
  if (!mng_ptr_->ms_one_ms_two_same_recal_) {
    return true;
  } else {
    // here we assume that precursor mass has same recal to MS2
    double prec_mass = deconv_ms_ptr_->getHeaderPtr()->getPrecMonoMass();
    double prec_ppo = (prec_mass  * (1 + recal) - refine_prec_mass_) / prec_mass;
    return std::abs(prec_ppo) <= ppo;
  }
}

// input is refineMsThree, result is score and recal and recalMass 
void ZeroPtmSlowMatch::compScore (ExtendMsPtr refine_ms_ptr, TheoPeakPtrVec theo_peaks,
                                  double ppo) {
  std::vector<double> ms_masses = getExtendMassVec(refine_ms_ptr);
  std::vector<double> theo_masses = getTheoMassVec(theo_peaks);
  std::vector<double> result_ppos = compMsMassPpos(ms_masses, theo_masses, ppo);

  if (!mng_ptr_->do_recalibration_) {
    recal_ = 0;
  } else {
    // minus is important 
    recal_ = -compAvg(result_ppos, mng_ptr_->recal_ppo_);
    if (!isValid(recal_, ppo)) {
      recal_ = 0;
    }
  }
  for (unsigned int i = 0; i < ms_masses.size(); i++) {
    ms_masses[i] = ms_masses[i] * (1 + recal_);
  }
  score_ = compNumMatchedTheoMasses(ms_masses, theo_masses, ppo);
}

// get result 
PrSMPtr ZeroPtmSlowMatch::geneResult() {
  SpParaPtr sp_para_ptr = mng_ptr_->prsm_para_ptr_->getSpParaPtr();
  return PrSMPtr(new PrSM(proteoform_ptr_, deconv_ms_ptr_, refine_prec_mass_, 
                          recal_, sp_para_ptr));
}

ZpSlowMatchPtrVec zeroPtmSlowFilter(DeconvMsPtr deconv_ms_ptr,
                                    ZpFastMatchPtrVec fast_matches,
                                    ZeroPtmMngPtr mng_ptr) {

  ZpSlowMatchPtrVec slow_matches;
  for (unsigned int i = 0; i < fast_matches.size(); i++) {
    ZpSlowMatchPtr slow_match = ZpSlowMatchPtr(
        new ZeroPtmSlowMatch(deconv_ms_ptr, fast_matches[i], mng_ptr));
    slow_matches.push_back(slow_match);
  }
  /* sort */
  std::sort(slow_matches.begin(), slow_matches.end(), compareZeroPtmSlowMatchDown);

  return slow_matches;
}

}
