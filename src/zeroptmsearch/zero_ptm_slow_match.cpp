#include "base/proteoform.hpp"
#include "base/activation.hpp"
#include "base/algorithm.hpp"
#include "spec/theo_peak.hpp"
#include "zeroptmsearch/zero_ptm_slow_match.hpp"

namespace prot {

ZeroPtmSlowMatch::ZeroPtmSlowMatch(int search_type,
                                   DeconvMsPtr deconv_ms_ptr, 
                                   ZpFastMatchPtr fast_match_ptr,
                                   ZeroPtmMngPtr mng_ptr) {
  mng_ptr_ = mng_ptr;
  deconv_ms_ptr_ = deconv_ms_ptr;
  fast_match_ptr_ = fast_match_ptr;

  proteoform_ptr_ = getSubProteoform(fast_match_ptr->getProteoformPtr(), 
                                     fast_match_ptr->getBegin(), 
                                     fast_match_ptr->getEnd());
  refine_prec_mass_ = proteoform_ptr_->getResSeqPtr()->getSeqMass();
  double delta = refine_prec_mass_ - deconv_ms_ptr->getHeaderPtr()->getPrecMonoMass();
  refine_ms_ptr_ = getMsThree(deconv_ms_ptr_, delta, mng_ptr_->sp_para_ptr_);

  ActivationPtr activation_ptr = deconv_ms_ptr_->getHeaderPtr()->getActivationPtr();
  NeutralLossPtr neu_loss_ptr = mng_ptr->base_data_ptr_->getNeutralLossNonePtr();
  double min_mass = mng_ptr_->sp_para_ptr_->getMinMass();
  TheoPeakPtrVec theo_peaks = getProteoformTheoPeak(proteoform_ptr_, 
                                                    activation_ptr, neu_loss_ptr, min_mass);

  compScore(refine_ms_ptr_, theo_peaks, mng_ptr_->sp_para_ptr_->getPeakTolerance()->getPpo());
}

// compute the average ppo
double compAvg(std::vector<double> ppos, double recal_ppo) {
  int cnt = 0;
  double sum = 0;
  for (unsigned int i = 0; i < ppos.size(); i++) {
    if (abs(ppos[i]) <= recal_ppo) {
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
    if (abs(prec_ppo) <= ppo) {
      return true;
    } else {
      return false;
    }
  }
}

// input is refineMsThree, result is score and recal and recalMass 
void ZeroPtmSlowMatch::compScore (ExtendMsPtr refine_ms_ptr, TheoPeakPtrVec theo_peaks,
                                  double ppo) {
  std::vector<double> ms_masses;
  getExtendMassVec(refine_ms_ptr, ms_masses);
  std::vector<double> theo_masses; 
  getTheoMassVec(theo_peaks, theo_masses);
  std::vector<double> result_ppos;
  compMsMassPpos(ms_masses, theo_masses, ppo, result_ppos);

  if (!mng_ptr_->do_recalibration_) {
    recal = 0;
  } else {
    // minus is important 
    recal = -compAvg(result_ppos, mng_ptr_->recal_ppo_);
    if (!isValid(recal, ppo)) {
      recal = 0;
    }
  }
  for (unsigned int i = 0; i < ms_masses.size(); i++) {
    ms_masses[i] = ms_masses[i] * (1 + recal);
  }
  score = compUniqueScore(ms_masses, theo_masses, ppo);
}

ZpSlowMatchPtrVec zeroPtmSlowFilter(int semi_align_type,
                                    DeconvMsPtr deconv_ms_ptr,
                                    ZpFastMatchPtrVec fast_matches,
                                    ZeroPtmMngPtr mng_ptr) {

  ZpSlowMatchPtrVec slow_matches;
  for (unsigned int i = 0; i < fast_matches.size(); i++) {
    ZpSlowMatchPtr slow_match = ZpSlowMatchPtr(
        new ZeroPtmSlowMatch(semi_align_type, deconv_ms_ptr, fast_matches[i], mng_ptr));
    slow_matches.push_back(slow_match);
  }
  /* sort */
  std::sort(slow_matches.begin(), slow_matches.end(), compareZeroPtmSlowMatchDown);

  return slow_matches;
}

}

/*
	public PrSM geneResult() throws Exception {
		// logger.debug(" gene result recal " + recal);
		return new PrSM(annoProtein, deconvMs, refinePrecMass, recal, mng.spPara);
	}


*/
