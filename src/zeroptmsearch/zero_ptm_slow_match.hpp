#ifndef ZERO_PTM_SLOW_MATCH_HPP_
#define ZERO_PTM_SLOW_MATCH_HPP_

#include "spec/deconv_peak.hpp"
#include "zeroptmsearch/zero_ptm_mng.hpp"
#include "zeroptmsearch/zero_ptm_fast_match.hpp"

namespace prot {

class ZeroPtmSlowMatch {
 public:
  ZeroPtmSlowMatch(DeconvMsPtr deconv_ms_ptr, ZpFastMatchPtr fast_match_ptr,
                   int search_type, ZeroPtmMngPtr mng_ptr);
	double getScore() {return score;}

 private:
	ZeroPtmMngPtr mng_ptr_;
  ZpFastMatchPtr fast_match_ptr_;
  DeconvMsPtr deconv_ms_ptr_;
  ProteoformPtr proteoform_ptr_;

	double refine_prec_mass_;
	ExtendMsPtr refine_ms_ptr_;

	double score = 0;
	double recal = 0;
};

}

#endif

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

	private void compScore(Ms<ExtendPeak> refineMsThree, TheoPeak ions[])
			throws Exception {
		double ppoDeviation[] = Pair.compPpoDeviation(refineMsThree, ions,
				mng.spPara.getPeakTolerance().getPpo());
		if (!mng.doRecalibration) {
			recal = 0;
		} else {
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
  */
