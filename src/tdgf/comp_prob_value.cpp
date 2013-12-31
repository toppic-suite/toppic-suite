#include <cmath>

#include "base/logger.hpp"
#include "tdgf/comp_prob_value.hpp"

namespace prot {

CompProbValue::CompProbValue(double convert_ratio, ResFreqPtrVec n_term_residues,
                             ResFreqPtrVec residues, int max_layer_num, 
                             int max_table_height, double max_sp_prec_mass) {
  convert_ratio_ = convert_ratio;
  n_term_acid_freq_sum_ = 0;
  for (unsigned int i = 0; i < n_term_residues.size(); i++) {
    int int_mass = (int)std::round(n_term_residues[i]->getMass() * convert_ratio);
    n_term_acid_masses_.push_back(int_mass);
    n_term_acid_frequencies_.push_back(n_term_residues[i]->getFreq());
    n_term_acid_freq_sum_ += n_term_residues[i]->getFreq();
  }
  double mass_sum = 0;
  double freq_sum = 0;
  for (unsigned int i = 0; i < residues.size(); i++) {
    int int_mass = (int)std::round(residues[i]->getMass() * convert_ratio);
    residue_masses_.push_back(int_mass);
    residue_frequencies_.push_back(residues[i]->getFreq());
    freq_sum = freq_sum + residues[i]->getFreq();
    mass_sum = mass_sum + residues[i]->getFreq() * residues[i]->getMass();
  }
  residue_avg_len_ = (int)std::round(mass_sum/freq_sum * convert_ratio);

  max_layer_num_ = max_layer_num;
  max_table_height_ = max_table_height;

  page_len_ = ORI_PAGE_LEN * (int) convert_ratio_;
  block_len_ = ORI_BLOCK_LEN * (int) convert_ratio_;
  page_table_ = new double[page_len_ * max_table_height_];

  max_sp_len_ = (int) std::round(max_sp_prec_mass * convert_ratio_);
  pos_scores_ = new short[max_sp_len_ + block_len_];
}

CompProbValue::~CompProbValue() {
  if (page_table_ != nullptr) {
    delete page_table_;
  }
  if (pos_scores_ != nullptr) {
    delete pos_scores_;
  }
}

void CompProbValue::compute(PrmPeakPtrVec peaks, int thresh, int shift_num, bool strict) {
  /*
  results = null; 
  // if score is less than 1, we do not compute p-value 
  if (thresh < 1) {
    return;
  }
  setMassErr(peaks, strict);
  setPosScores(peakMasses, peakTolerances, baseTypes);
  int lastPeakIndex = peakMasses.length - 1;
  int maxPeakMass = peakMasses[lastPeakIndex]
      + peakTolerances[lastPeakIndex];
  setHeight(thresh, maxPeakMass);
  setPeakBgnEnd(peakMasses, peakTolerances);
  if (spLen + residueAvgLen >= maxSpLen) {
    logger.error("Spectral precursor mass is too large"
                 + peaks[peaks.length - 1].getMonoMass());
    return;
  }
  if (nShift + 1 > maxLayerNum) {
    logger.error("Number of unexpected PTMs is too large" + nShift);
    return;
  }
  this.nShift = nShift;
  comp();
  */
}

void CompProbValue::setMassErr(PrmPeakPtrVec peaks, bool strict) {
  LOG_DEBUG("MAX MASS " << peaks[peaks.size() - 1]->getMonoMass());
  for (unsigned int i = 0; i < peaks.size(); i++) {
    peak_masses_.push_back((int) std::round(peaks[i]->getMonoMass()
                                            * convert_ratio_));
    // we use NStrict and CRelax tolerance
    if (strict) {
      peak_tolerances_.push_back(std::ceil(peaks[i]->getStrictTolerance()
                                           * convert_ratio_));
    } else {
      peak_tolerances_.push_back(std::ceil(peaks[i]->getNStrictCRelaxTolerance()
                                           * convert_ratio_));
    }
    base_types_.push_back(peaks[i]->getBaseType());
  }
  // tolerance is 0 for mass 0 
  peak_tolerances_[0] = 0;
  // we use restrict error tolerance for residue masses 
  int last_peak_index = peaks.size() - 1;
  peak_tolerances_[last_peak_index] = std::floor(peaks[last_peak_index]->getStrictTolerance()
                                                * convert_ratio_);
}

void CompProbValue::setPosScores(std::vector<int> &peak_masses, 
                                 std::vector<int> &peak_tolerances,
                                 std::vector<int> &base_types) {
  int len = sizeof(pos_scores_)/sizeof(short); 
  LOG_DEBUG("pos scr length " << len);
  memset(pos_scores_, 0, sizeof(pos_scores_));
  // mass 0 and residue sum mass are not used for scoring 
  for (unsigned int i = 1; i < peak_masses.size() - 1; i++) {
    // here we use ceil/floor not round since each unit is very accurate
    int bgn = peak_masses[i] - peak_tolerances[i];
    if (bgn < 0) {
      bgn = 0;
    }
    int end = peak_masses[i] + peak_tolerances[i];
    if (end >= len) {
      end = len - 1;
    }
    if (base_types[i] == PRM_PEAK_TYPE_ORIGINAL) {
      // the peak is m
      for (int p = bgn; p <= end; p++) {
        if (pos_scores_[p] == 0) {
          pos_scores_[p] = 1;
        }
        if (pos_scores_[p] == 2) {
          pos_scores_[p] = 3;
        }
      }
    } else {
      for (int p = bgn; p <= end; p++) {
        if (pos_scores_[p] == 0) {
          pos_scores_[p] = 2;
        }
        if (pos_scores_[p] == 1) {
          pos_scores_[p] = 3;
        }
      }
    }

  }
  for (int i = 0; i < len; i++) {
    if (pos_scores_[i] == 1 || pos_scores_[i] == 2) {
      pos_scores_[i] = 1;
    } else if (pos_scores_[i] == 3) {
      pos_scores_[i] = 2;
    }
  }
}

void CompProbValue::setHeight(int thresh, int max_peak_mass) {
  height_ = thresh + 1;
  if (height_ > max_table_height_) {
    height_ = max_table_height_;
  }
  sp_len_ = max_peak_mass + residue_avg_len_ + 2;
  sp_table_size_ = sp_len_ * height_;
  page_table_size_ = page_len_ * height_;
  block_table_size_ = block_len_ * height_;
  LOG_DEBUG("max mass " << max_peak_mass << " sp len " << sp_len_
            << " blk size " << block_table_size_ << " page size " << page_table_size_);
  LOG_DEBUG("height " << height_);
  acid_dists_.clear();
  for (unsigned int i = 0; i < residue_masses_.size(); i++) {
    acid_dists_.push_back(residue_masses_[i] * height_);
  }
}

}

/*
	 void setPeakBgnEnd(int peakMasses[], int peakTolerances[]) {
		peakMassBgns = new int[peakMasses.length];
		peakMassEnds = new int[peakMasses.length];
		peakTableBgns = new int[peakMasses.length];
		peakTableEnds = new int[peakMasses.length];
		for (int i = 0; i < peakMasses.length; i++) {
			peakMassBgns[i] = peakMasses[i] - peakTolerances[i];
			peakMassEnds[i] = peakMasses[i] + peakTolerances[i];
			peakTableBgns[i] = peakMassBgns[i] * height;
			peakTableEnds[i] = peakMassEnds[i] * height + (height - 1);
		}
	}

	 void comp() {
		results = new double[nShift + 1][][];
		priors = new double[nShift + 1][][];
		double oneLayerResults[][][];
        logger.debug("Start computation");
		oneLayerResults = compOneLayer(null, null, true);
		// one layer results [0] is score probability for each peak
		results[0] = oneLayerResults[0]; 
		// one layer results [1] is prior probability for next layer
		priors[0] = oneLayerResults[1]; 
		for (int i = 1; i < nShift + 1; i++) {
			oneLayerResults = compOneLayer(results[i - 1], priors[i - 1], false);
			results[i] = oneLayerResults[0];
			priors[i] = oneLayerResults[1];
		}
		factors = compFactors();
	}
	
	 double[] compFactors() {
		double factors[] = new double[nShift + 1];
		// zero ptm
		int lastPeakIndex = peakMasses.length - 1;
		double hitSum = 0;
		for (int i = 0; i < height; i++) {
			hitSum += results[0][lastPeakIndex][i];
		}
		if (hitSum == 0) {
			hitSum = 1;
		}
		factors[0] = 1 /hitSum;

//		if (nShift >= 1) {
//			double priorSum = 0;
//			for (int i = 0; i < height; i++) {
//				priorSum += priors[0][lastPeakIndex][i];
//			}
//			// one ptm
//			factors[1] = 1/priorSum;
//			// divided by the last peakEnd - last peakBegin
//			double peakWidth = peakMassEnds[lastPeakIndex] - peakMassBgns[lastPeakIndex] + 1;
//			factors[1] = factors[1] / peakWidth;
//			factors[1] = factors[1] * K;
//			// two ptm
//			for (int i = 2; i < nShift+1; i++) {
//				factors[i] = 1 / (priorSum * i);
//				// the reason is that peakWidth is used once for each mass shift 
//				factors[i] = factors[i] / Math.pow(peakWidth, i);
//				if (i == 2) {
//					factors[i] = factors[i] * K2;
//				}
//				else {
//					factors[i] = factors[i] * Math.pow(K,i);
//				}
//			}
//		}
		
		double priorSum = 0;
		for (int i = 0; i < height; i++) {
			priorSum += priors[0][lastPeakIndex][i];
		}
		double peakWidth = peakMassEnds[lastPeakIndex] - peakMassBgns[lastPeakIndex] + 1;
		if (nShift == 1) {
			// one ptm
			factors[1] = 1/priorSum;
			// divided by the last peakEnd - last peakBegin
			factors[1] = factors[1] / peakWidth;
			factors[1] = factors[1] * K;
		}
		if (nShift == 2) {
			// two ptm
			factors[2] = 1/ (priorSum * 2);
			// the reason is that peakWidth is used once for each mass shift 
			factors[2] = factors[2] / Math.pow(peakWidth, 2);
			factors[2] = factors[2] * K2;
		}
		return factors;
	}

	 double[][][] compOneLayer(double[][] results, double[][] priors,
			boolean isFirstLayer) {
		double probs[][][] = new double[2][peakMasses.length][height];
		Arrays.fill(pageTable, 0);
		int peakIndex = 0;
		for (int winTableBgn = 0; winTableBgn < spTableSize; winTableBgn = winTableBgn
				+ blockTableSize) {
			// clear result block
			int WinTableEnd = winTableBgn + blockTableSize - 1;
			int pagePos = winTableBgn % pageTableSize;
			runClear(pagePos);
			if (isFirstLayer) {
				runFirstLayerInit(winTableBgn, WinTableEnd);
			} else {
				runInit(results, priors, winTableBgn, WinTableEnd, peakIndex);
			}
			for (int i = 0; i < acidDists.length; i++) {
				if (residueFrequencies[i] <= 0) {
					continue;
				}
				int prevBlockBgn = pagePos - acidDists[i];
				if (prevBlockBgn < 0) {
					prevBlockBgn += pageTableSize;
				}
				int prevBlockEnd = prevBlockBgn + blockTableSize;
				if (prevBlockEnd < pageTableSize) {
					runAddProb(pagePos, prevBlockBgn, blockTableSize,
							residueFrequencies[i]);
				} else {
					int firstPartSize = pageTableSize - prevBlockBgn;
					runAddProb(pagePos, prevBlockBgn, firstPartSize,
							residueFrequencies[i]);
					runAddProb(pagePos + firstPartSize, 0, blockTableSize
							- firstPartSize, residueFrequencies[i]);
				}
			}
			int spBgn = winTableBgn / height;
			int spEnd = winTableBgn / height + blockLen - 1;
			runUpdate(pagePos + blockTableSize - 1, spBgn, spEnd);
			// update peakPnt and get probs 
			while (peakIndex < peakMasses.length
					&& peakTableEnds[peakIndex] <= WinTableEnd) {
				for (int i = peakMassBgns[peakIndex]; i <= peakMassEnds[peakIndex]; i++) {
					for (int j = 0; j < height; j++) {
						int pos = (i * height + j) % pageTableSize;
						probs[0][peakIndex][j] += pageTable[pos];
					}
				}

				int priorMassBgn = peakMassBgns[peakIndex] - residueAvgLen;
				int priorMassEnd = peakMassBgns[peakIndex] - 1;
				// double sum = 0;
				for (int i = priorMassBgn; i <= priorMassEnd; i++) {
					for (int j = 0; j < height; j++) {
						int pos = (i * height + j) % pageTableSize;
						if (pos < 0) {
							pos = pos + pageTableSize;
						}
						probs[1][peakIndex][j] += pageTable[pos];
					}
				}
				peakIndex++;
			}
		}
		return probs;
	}

	 void runClear(int pagePos) {
		for (int i = pagePos; i < pagePos + blockTableSize; i++) {
			pageTable[i] = 0f;
		}
	}

	 void runFirstLayerInit(int winTableBgn, int winTableEnd) {
		double baseProb = 1.0;
		for (int i = 0; i < nTermAcidMasses.length; i++) {
			if (nTermAcidFrequences[i] <= 0) {
				continue;
			}
			//logger.debug("nTermAcidMass " + nTermAcidMasses[i] + " freq " +
			// nTermAcidFrequences[i]);
			int pos = nTermAcidMasses[i];
			int k = pos * height;
			// logger.debug("k " + k + " end " + blockTableSize);
			if (k >= winTableBgn && k <= winTableEnd) {
				pageTable[k % pageTableSize] += baseProb
						* nTermAcidFrequences[i];
				// logger.debug("pageTable value " + pageTable[k %
				// pageTableSize]);
			}
		}
	}

	 void runInit(double results[][], double priors[][],
			int winTableBgn, int winTableEnd, int peakIndex) {
		// zero
		double baseProb = 1.0;
		for (int i = 0; i < nTermAcidMasses.length; i++) {
			if (nTermAcidFrequences[i] <= 0) {
				continue;
			}
			int pos = nTermAcidMasses[i];
			for (int k = pos * height; k < pos * height + residueAvgLen
					* height; k += height) {
				if (k >= winTableBgn && k <= winTableEnd) {
					pageTable[k % pageTableSize] += baseProb
							* nTermAcidFrequences[i];
				}
				if (k > winTableEnd) {
					break;
				}
			}
		}

		// results
		for (int i = 0; i < peakMasses.length; i++) {
			if (peakTableEnds[i] >= winTableEnd) {
				break;
			}
			if (peakTableEnds[i] + residueAvgLen * height > winTableBgn) {
				for (int k = peakTableEnds[i] + 1; k < peakTableEnds[i]
						+ residueAvgLen * height; k += height) {
					if (k >= winTableBgn && k <= winTableEnd) {
						for (int h = 0; h < height; h++) {
							pageTable[(k + h) % pageTableSize] += results[i][h];
						}
					}
					if (k > winTableEnd) {
						break;
					}
				}
			}
		}
	}

	 void runAddProb(int pagePos, int prevBlockBgn, int size, double f) {
		int prev_pos = prevBlockBgn;
		for (int i = pagePos; i < pagePos + size; i++) {
			pageTable[i] = pageTable[i] + pageTable[prev_pos] * f;
			prev_pos++;
		}
	}

	 void updateCol(int colEnd, int scr) {
		int colBgn = colEnd - height + 1;
		int j, pre;
		for (j = 1; j <= scr; j++) {
			pre = colEnd - j;
			if (pre >= colBgn) {
				pageTable[colEnd] += pageTable[pre];
			}
		}
		for (j = colEnd - scr - 1; j >= colBgn; j--) {
			pageTable[j + scr] = pageTable[j];
		}
		for (j = colBgn; j < colBgn + scr; j++) {
			if (j <= colEnd) {
				pageTable[j] = 0;
			}
		}
	}

	 void runUpdate(int pageEnd, int scrBgn, int scrEnd) {
		int colEnd = pageEnd;
		logger.trace("scr bgn " + scrBgn + " end " + scrEnd);
		for (int i = scrEnd; i >= scrBgn; i--) {
			int scr = posScores[i];
			if (scr > 0) {
				updateCol(colEnd, scr);
			}
			colEnd = colEnd - height;
		}
	}

	
	 double getRawProb(int shift, int thresh) {
		double probSum = 0;
		if (results == null) {
			return 1.0;
		}
        int lastPeakIndex = peakMasses.length - 1;
        for (int score = thresh; score < height; score++) {
            probSum = probSum + results[shift][lastPeakIndex][score];
        }
        //System.out.println("Shift " + shift + " thresh " + thresh + " probSum " + probSum + " factor " + factors[shift] + " prob " + (probSum * factors[shift]));
        return probSum * factors[shift];
	}
	
	// The difference between getProb and getRawProb is nTermAcidFreqSum
	 double getProb(int shift, int thresh) {
        System.out.println("get prob");
		if (results == null) {
			return 1;
		}
		double rawProb = getRawProb(shift, thresh);
		return rawProb * nTermAcidFreqSum;
	}

	 double getOneValueProb(int shift, int value) {
		if (results == null) {
			return 1.0;
		}
		int lastPeakIndex = peakMasses.length - 1;
		double prob = results[shift][lastPeakIndex][value] * factors[shift];
		return prob * nTermAcidFreqSum;
	}
}
*/
