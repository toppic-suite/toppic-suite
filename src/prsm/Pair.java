package edu.ucsd.proteomics.prsm.pair;

import java.util.ArrayList;
import java.util.Arrays;

import org.apache.log4j.Logger;

import com.jap.proteomics.base.ion.Ion;
import com.jap.proteomics.base.util.BioArray;
import com.jap.proteomics.spec.extendsp.ExtendPeak;
import com.jap.proteomics.spec.peak.PositionComparator;
import com.jap.proteomics.spec.sp.Ms;
import com.jap.proteomics.spec.theosp.TheoPeak;



public class Pair implements Comparable<Pair> {
    
    public static Logger logger = Logger.getLogger(Pair.class);

	protected int x, y;

	public Pair(int x, int y) {
		this.x = x;
		this.y = y;
	}

	public int getX() {
		return x;
	}

	public int getY() {
		return y;
	}

	public void setX(int x) {
		this.x = x;
	}

	public void setY(int y) {
		this.y = y;
	}

	/** compare the pairs based on x and y */
	public int compareTo(Pair p2) {
		int x2 = p2.getX();
		int y2 = p2.getY();
		if (y2 - y < 0) {
			return 1;
		} else if (y2 - y > 0) {
			return -1;
		} else if (x2 - x < 0) {
			return 1;
		} else if (x2 - x > 0) {
			return -1;
		} else {
			return 0;
		}
	}

	public static boolean increaseIJ(int i, int j, double deviation, 
			double tolerance, double msMasses[], double seqMasses[]) {
		/*
		 * we assume that each real peak is matched to at most one theoretical
		 * peak, so we do not check i and j+1
		 */
		if (deviation <= 0) {
			return true;
		}
		/* severl real peak can be matched to the same theoretical peak */
		if (i >= msMasses.length - 1) {
			return false;
		}

		double nextPos = msMasses[i+1];
		if (Math.abs(nextPos - seqMasses[j]) <=  tolerance
				&& (j == seqMasses.length - 1 || Math.abs(nextPos
						- seqMasses[j]) < Math.abs(nextPos - seqMasses[j + 1]))) {
			return true;
		} else {
			return false;
		}
	}
	
	/** compute deviation for each peak */
	private static double[] compPpoDeviation(double msMasses[],
	            double theoMasses[], double ppo) throws Exception {
	        double minDistances[] = new double[msMasses.length];
	        Arrays.fill(minDistances, Double.MAX_VALUE);
	        /* extendMsThree do not have 0 and precursor mass */
	        int i = 0;
	        int j = 0;
	        while (i < msMasses.length && j < theoMasses.length) {
	            double d = msMasses[i] - theoMasses[j];
	            if (Math.abs(d) <= Math.abs(minDistances[i])) {
	                minDistances[i] = d;
	            }
	            double tolerance = msMasses[i] * ppo;
	            if (increaseIJ(i, j, d,  tolerance, msMasses, theoMasses)) {
	                i++;
	            } else {
	                j++;
	            }
	        }
	        // change distance to ppo
	        double minPpos[] = new double[msMasses.length];
	        for (i = 0; i < msMasses.length; i++) {
	            if (msMasses[i] > 0) {
	                minPpos[i] = minDistances[i] / msMasses[i];
	            }
	            else {
	                minPpos[i] = Double.MAX_VALUE;
	            }
	        }
	        return minPpos;
	    }
	
	
	/** compute unique score */
    private static double compUniqScore(double msMasses[],
            double theoMasses[], double ppo) throws Exception {
        double minDistances[] = new double[theoMasses.length];
        Arrays.fill(minDistances, Double.MAX_VALUE);
        /* extendMsThree do not have 0 and precursor mass */
        int i = 0;
        int j = 0;
        while (i < msMasses.length && j < theoMasses.length) {
            double d = msMasses[i] - theoMasses[j];
            if (Math.abs(d) <= Math.abs(minDistances[j])) {
                minDistances[j] = d;
            }
            double tolerance = msMasses[i] * ppo;
            if (increaseIJ(i, j, d, tolerance, msMasses, theoMasses)) {
                i++;
            } else {
                j++;
            }
        }
        double score = 0;
        for (i = 0; i < theoMasses.length; i++) {
            if (theoMasses[i] > 0) {
                double error = minDistances[i] / theoMasses[i];
                if (Math.abs(error) <= ppo) {
                    score += 1.0;
                }
            }
        }
        return score;
    }
	
    /******************************************************************
     * deviation computation for each experimental peak. used in 
     * zero_ptm PrSM recalibration. 
     ******************************************************************/
	public static double[] compPpoDeviation(Ms<ExtendPeak> ms,
			TheoPeak ions[], double ppo) throws Exception {
		double ionMasses[] = new double[ions.length];
		for (int k = 0; k < ions.length; k++) {
			ionMasses[k] = ions[k].getModMass();
		}
		double msMasses[] = ms.getPositions();
		return compPpoDeviation(msMasses, ionMasses, ppo);
	}
	
	/******************************************************************
	 * Diagonal score computation used in zero_ptm PrSM, 
	 * MinDistance is for ions, not experimental peaks
	 ******************************************************************/
	public static double compIonScore(Ms<ExtendPeak> ms,
			TheoPeak ions[], double recal, double ppo) throws Exception {
		double ionMasses[] = new double[ions.length];
		for (int k = 0; k < ions.length; k++) {
			ionMasses[k] = ions[k].getModMass();
		}
		double msMasses[] = ms.getPositions();
		for (int i = 0; i < msMasses.length; i++) {
		    msMasses[i] *= (1 + recal);
		}
		return compUniqScore(msMasses, ionMasses, ppo);
	}
	
    /***************************************************************
     * AnnoProtein factory
     * *************************************************************/
    /** get all the matching points on the best diagonal */
    public static PeakIonPair[] findPairs(Ms<ExtendPeak> alignMs,
            TheoPeak theoPeaks[], int bgn, int end) throws Exception {

        ArrayList<PeakIonPair> pairList = new ArrayList<PeakIonPair>();
        Arrays.sort(theoPeaks, new PositionComparator());
        double ionMasses[] = new double[theoPeaks.length];
        for (int k = 0; k < theoPeaks.length; k++) {
            ionMasses[k] = theoPeaks[k].getModMass();
        }
        double realMasses[] = alignMs.getPositions();
        int i = 0;
        int j = 0;
        while (i < alignMs.size() && j < theoPeaks.length) {

            double deviation = alignMs.getPosition(i) - theoPeaks[j].getModMass();
            Ion ion = theoPeaks[j].getIon();
            double err = alignMs.get(i).getOrigTolerance();
            //logger.debug("i " + i  + " " + alignMs.getPosition(i) + " j " + theoPeaks[j].getIon().getDisplayName() + " " + theoPeaks[j].getModMass() + " " + deviation + " err " + err);
            if (ion.getPos() >= bgn && ion.getPos() <= end) {
                if (Math.abs(deviation) <= err) {
                    //logger.debug("added");
                    pairList.add(new PeakIonPair(alignMs.get(i), theoPeaks[j]));
                }
            }
            if (increaseIJ(i, j, deviation, err, realMasses, ionMasses)) {
                i++;
            } else {
                j++;
            }
        }
        return (PeakIonPair[]) pairList.toArray(new PeakIonPair[0]);
    }
    
    /***************************************************************
     * Precursor mass correction, called by PtmSlowMatch
     * *************************************************************/
    public static double[] getNCScore(Ms<ExtendPeak> alignMs,
            TheoPeak theoPeaks[], int bgn, int end, double delta, double ppo) throws Exception {
        double msMasses[] = alignMs.getPositions();
        
        // n score
        ArrayList<Double> theoNMassList = new ArrayList<Double> ();
        for (int k = 0; k < theoPeaks.length; k++) {
            Ion ion = theoPeaks[k].getIon();
            int pos = ion.getPos();
            if (ion.getIonType().isNTerm() && pos >= bgn && pos <= end ) {
                theoNMassList.add(theoPeaks[k].getModMass() + delta);
            }
        }
        double theoMasses[] = BioArray.cnvtDoubleList(theoNMassList);
        double nScore = compUniqScore(msMasses, theoMasses, ppo);
        
        // c score 
        ArrayList<Double> theoCMassList = new ArrayList<Double> ();
        for (int k = 0; k < theoPeaks.length; k++) {
            Ion ion = theoPeaks[k].getIon();
            int pos = ion.getPos();
            //logger.debug("is n term " + ion.getIonType().isNTerm() + " pos " + pos + " mass " + theoPeaks[k].getModMass());
            if (!ion.getIonType().isNTerm() && pos >= bgn && pos <= end ) {
                theoCMassList.add(theoPeaks[k].getModMass() + delta);
            }
        }
        theoMasses = BioArray.cnvtDoubleList(theoCMassList);
        double cScore = compUniqScore(msMasses, theoMasses, ppo);
        
        double result[] = new double[2];
        result[0] = nScore;
        result[1] = cScore;
        return result;
    }
}
