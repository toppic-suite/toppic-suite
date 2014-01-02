#include "tdgf/count_test_num.hpp"

namespace prot {

CountTestNum::CountTestNum() {
}

}
/*
	public CountTestNum(BpSpec seqs[], NModBpSpec nModSeqs[], ResFreq nTermResidues[], ResFreq residues[], TdgfMng mng) {
		this.mng = mng;
		this.seqs = seqs;
		this.nModSeqs = nModSeqs;
		convertRatio = mng.doubleToIntConstant;
		maxSpLen = (int) Math.round(mng.maxSpPrecMass * convertRatio);
		//System.out.println("max prec mass " + mng.maxSpPrecMass + " max sp len " + maxSpLen);
		residueAvgLen = (int)(ResFreqArrayUtil.getAvgMass(residues) * convertRatio);
		normFactor = ArrayUtil.getSum(ResFreqArrayUtil.getFrequencies(nTermResidues));
		initCompMassCnt(nModSeqs);
		initPrefMassCnt(nModSeqs);
		initSuffMassCnt(seqs);
		initInternalMassCnt();
	}
  */

	/** initialize the four tables for mass counts */
/*
	private void initCompMassCnt(NModBpSpec nModSeqs[]) {
		nonPtmCompMassCnts = new double[maxSpLen]; 
		for (int i = 0; i < nModSeqs.length; i++) {
			double m = nModSeqs[i].getResSeq().getResMassSum();
			//System.out.println("mass " + m);
			nonPtmCompMassCnts[convertMass(m)] += 1.0;
		}
	}
*/
	/** initialize the four tables for mass counts */
/*
	private void initPrefMassCnt(NModBpSpec nModSeqs[]) {
		nonPtmPrefMassCnts = new double[maxSpLen];
		for (int i = 0; i < nModSeqs.length; i++) {
			double extBMasses[] = nModSeqs[i].getExtBMasses();
			// prefix
			for (int j = 1; j < extBMasses.length - 1; j++) {
				nonPtmPrefMassCnts[convertMass(extBMasses[j])] += 1.0;
			}
		}
	}
  */

	/** initialize the four tables for mass counts */
/*
	private void initSuffMassCnt(BpSpec seqs[]) {
		// sequence mass 
		nonPtmSuffMassCnts = new double[maxSpLen];
		for (int i = 0; i < seqs.length; i++) {
			BreakPoint[] extBps = seqs[i].getExtBps();
			// suffix
			for (int j = 1; j < extBps.length - 1; j++) {
				nonPtmSuffMassCnts[convertMass(extBps[j].getSrm())] += 1;
			}
		}
		//nonPtmSuffMassCnts = ptmSuffMassCnts;
	}

	private void initInternalMassCnt() {
		nonPtmInternalMassCnts = new double[maxSpLen];
		// middle
		double normCount = 0;
		// use approxiation to speed up
		for (int i = maxSpLen - 1; i >= 0; i--) {
			normCount += nonPtmSuffMassCnts[i];
			nonPtmInternalMassCnts[i] = normCount / residueAvgLen;
		}
		//nonPtmInternalMassCnts = ptmInternalMassCnts;
	}

	private int convertMass(double m) {
		int n = (int) Math.round(m * convertRatio);
		if (n < 0) {
			logger.warn("Negative mass value : " + m);
			return 0;
		}
		if (n >= maxSpLen) {
			n = maxSpLen - 1;
		}
		return n;
	}

	public double compCandNum(SemiAlignType type, int nShift, double oriMass,
			double oriTolerance) {
		double nCandidates = 0;
		if (nShift == 0) {
			nCandidates = compNormNonPtmCandNum(type, nShift, oriMass, oriTolerance);
		}
		// with shifts 
		else if (nShift == 1){
			nCandidates = compOnePtmCandNum(type, nShift, oriMass);
		}
		else {
			nCandidates = compMultiplePtmCandNum(type, nShift, oriMass);
		}

		if (nCandidates == 0.0) {
			logger.warn("candidate number is ZERO");
		}
		if (type == SemiAlignType.PREFIX || type == SemiAlignType.SUFFIX) {
			nCandidates = nCandidates * PREFIX_SUFFIX_ADJUST;
		}
		else if (type == SemiAlignType.INTERNAL) {
			nCandidates = nCandidates * INTERNAL_ADJUST;
		}
    	return nCandidates;
	}
	
	private double compNormNonPtmCandNum(SemiAlignType type, int nshift, double oriMass, double oriTolerance) {
		int low = (int) Math.floor((oriMass - oriTolerance) * convertRatio);
		int high = (int) Math.ceil((oriMass + oriTolerance) * convertRatio);
		double nCandidates = compSeqNum(type, low, high);
		// normalization: the reason is that we a residue list with frequency sum > 1 in CompProbValue 
		if (type == SemiAlignType.COMPLETE || type == SemiAlignType.PREFIX) {
			nCandidates = nCandidates /normFactor;
			//System.out.println("nCandidate " + nCandidates + " normFactor " + normFactor);
		} 
		return nCandidates;
	}
	
	private double compOnePtmCandNum (SemiAlignType type, int nShift, double oriMass ) {
		double nCandidates = 0;
		if (type == SemiAlignType.COMPLETE) {
			nCandidates = seqs.length;
		} else if (type == SemiAlignType.PREFIX) {
			for (int i = 0; i < seqs.length; i++) {
				nCandidates += seqs[i].getResSeq().getLen();
			}
		} else if (type == SemiAlignType.SUFFIX) {
			for (int i = 0; i < seqs.length; i++) {
				nCandidates += seqs[i].getResSeq().getLen();
			}
		} else if (type == SemiAlignType.INTERNAL) {
			for (int i = 0; i < seqs.length; i++) {
				nCandidates = nCandidates + seqs[i].getResSeq().getLen() * seqs[i].getResSeq().getLen();
			}
		}
		return nCandidates;
	}
	
	private double compMultiplePtmCandNum (SemiAlignType type, int nShift, double oriMass ) {
		double nCandidates = 0;
		if (type == SemiAlignType.COMPLETE) {
			for (int i = 0; i < seqs.length; i++) {
				nCandidates = nCandidates
						+ Math.pow(seqs[i].getResSeq().getLen(), nShift - 1);
			}
		} else if (type == SemiAlignType.PREFIX) {
			for (int i = 0; i < seqs.length; i++) {
				nCandidates += Math.pow(seqs[i].getResSeq().getLen(), nShift);
			}
		} else if (type == SemiAlignType.SUFFIX) {
			for (int i = 0; i < seqs.length; i++) {
				nCandidates += Math.pow(seqs[i].getResSeq().getLen(), nShift);
			}
		} else if (type == SemiAlignType.INTERNAL) {
			for (int i = 0; i < seqs.length; i++) {
				nCandidates = nCandidates
						+ Math.pow(seqs[i].getResSeq().getLen(), nShift + 1);
			}
		}
		// use average proportion to estimate the number of matched substring. 
		nCandidates = nCandidates * getAvgProportion(oriMass, mng.spPara.getPeakTolerance().getPpo(), convertRatio, residueAvgLen);
		return nCandidates;
	}

	private double compSeqNum(SemiAlignType type, int low, int high) {
		double candNum = 0;
		if (type == SemiAlignType.COMPLETE) {
			candNum = compMassNum(nonPtmCompMassCnts, low, high);
		} else if (type == SemiAlignType.PREFIX) {
			candNum = compMassNum(nonPtmPrefMassCnts, low, high);
		} else if (type == SemiAlignType.SUFFIX) {
			candNum = compMassNum(nonPtmSuffMassCnts, low, high);
		} else if (type == SemiAlignType.INTERNAL) {
			candNum = compMassNum(nonPtmInternalMassCnts, low, high);
		}
		return candNum;
	}

	private double compMassNum(double cnts[], int low, int high) {
		double cnt = 0;
		if (high >= maxSpLen) {
			high = maxSpLen - 1;
		}
		if (low < 0) {
			low = 0;
		}
		if (low > high) {
			low = high;
		}
		for (int i = low; i <= high; i++) {
			cnt += cnts[i];
		}
		return cnt;
	}
	
	public static double getAvgProportion(double mass, double ppo, double convertRatio, double residueAvgLen ) {
		double proportion = mass * ppo * 2  * convertRatio/ residueAvgLen;
		if (proportion == 0) {
			logger.error("Error in computing proportion = 0");
		}
		return proportion;
	}
}
*/
