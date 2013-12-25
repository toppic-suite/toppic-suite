#include "base/proteoform.hpp"
#include "base/activation.hpp"
#include "spec/theo_peak.hpp"
#include "zeroptmsearch/zero_ptm_slow_match.hpp"

namespace prot {

ZeroPtmSlowMatch::ZeroPtmSlowMatch(DeconvMsPtr deconv_ms_ptr, 
                                   ZpFastMatchPtr fast_match_ptr,
                                   int search_type, 
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

  /*
  compScore(refineMsThree, theoPeaks);
  */
}

}

/*
	public PrSM geneResult() throws Exception {
		// logger.debug(" gene result recal " + recal);
		return new PrSM(annoProtein, deconvMs, refinePrecMass, recal, mng.spPara);
	}

	public int compareTo(ZeroPtmSlowMatch y) {
		double y_scr = y.getScore();
		if (y_scr - getScore() < 0) {
			return -1;
		} else if (y_scr - getScore() > 0) {
			return 1;
		} else {
			return 0;
		}
	}

	// private functions for recalibration

  // input is refineMsThree, result is score and recal and recalMass 
	private void compScore(Ms<ExtendPeak> refineMsThree, TheoPeak ions[])
			throws Exception {
		double ppoDeviation[] = Pair.compPpoDeviation(refineMsThree, ions,
				mng.spPara.getPeakTolerance().getPpo());
		if (!mng.doRecalibration) {
			recal = 0;
		} else {
			// minus is important 
			recal = -compAvg(ppoDeviation, mng.recalPpo);
			if (!isValid(recal)) {
				recal = 0;
			}
		}
		score = Pair.compIonScore(refineMsThree, ions, recal, mng.spPara.getPeakTolerance().getPpo());
		// logger.debug("Recalibration " + recal + " score " + score);
	}

	private double compAvg(double d[], double recal_shift_width) {
		int cnt = 0;
		double sum = 0;
		for (int i = 0; i < d.length; i++) {
			if (Math.abs(d[i]) <= recal_shift_width) {
				cnt++;
				sum += d[i];
			}
		}
		if (cnt == 0) {
			return 0;
		} else {
			return sum / cnt;
		}
	}
*/
	/**
	 * compute the validation of the candidates based on the different of
	 * precursor mass
	 */

/*
	private boolean isValid(double recal) {
		if (!mng.isMs1Ms2SameRecal) {
			return true;
		} else {
			// here we assume that precursor mass has same recal to MS2
			double spMass = deconvMs.getHeader().getPrecMonoMass();
			double ppo = (spMass * (1 + recal) - refinePrecMass) / spMass;
			if (Math.abs(ppo) <= mng.spPara.getPeakTolerance().getPpo()) {
				return true;
			} else {
				return false;
			}
		}
	}
}
*/
