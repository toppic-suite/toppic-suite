#include <cmath>
#include <cstring>

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
  results_.clear(); 
  // if score is less than 1, we do not compute p-value 
  if (thresh < 1) {
    return;
  }
  setMassErr(peaks, strict);
  setPosScores(peak_masses_, peak_tolerances_, base_types_);
  int last_peak_index = peak_masses_.size() - 1;
  int max_peak_mass = peak_masses_[last_peak_index]
      + peak_tolerances_[last_peak_index];
  setHeight(thresh, max_peak_mass);
  setPeakBgnEnd(peak_masses_, peak_tolerances_);
  if (sp_len_ + residue_avg_len_ >= max_sp_len_) {
    LOG_ERROR("Spectral precursor mass is too large "
              << peaks[peaks.size() - 1]->getMonoMass());
    return;
  }
  if (shift_num + 1 > max_layer_num_) {
    LOG_ERROR("Number of unexpected PTMs is too large " << shift_num);
    return;
  }
  shift_num_ = shift_num;
  comp();
}

void CompProbValue::setMassErr(PrmPeakPtrVec &peaks, bool strict) {
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

void CompProbValue::setPeakBgnEnd(std::vector<int> &peak_masses, 
                                  std::vector<int> &peak_tolerances) {
  peak_mass_bgns_.clear();
  peak_mass_ends_.clear();
  peak_table_bgns_.clear();
  peak_table_ends_.clear();
  for (unsigned int i = 0; i < peak_masses.size(); i++) {
    peak_mass_bgns_.push_back(peak_masses[i] - peak_tolerances[i]);
    peak_mass_ends_.push_back(peak_masses[i] + peak_tolerances[i]);
    peak_table_bgns_.push_back(peak_mass_bgns_[i] * height_);
    peak_table_ends_.push_back(peak_mass_ends_[i] * height_ + (height_ -1));
  }
}

void CompProbValue::comp() {
  results_.clear();
  priors_.clear();
  std::vector<std::vector<double>> one_layer_results;
  std::vector<std::vector<double>> one_layer_priors;
  LOG_DEBUG("Start computation");
  std::vector<std::vector<double>> empty_vec;
  compOneLayer(empty_vec, true, one_layer_results, one_layer_priors);
  // score probability for each peak
  results_.push_back(one_layer_results);
  // prior probability for next layer
  priors_.push_back(one_layer_priors);
  for (int i = 1; i < shift_num_ + 1; i++) {
    one_layer_results.clear();
    one_layer_priors.clear();
    compOneLayer(results_[i - 1], false, one_layer_results, one_layer_priors);
    results_.push_back(one_layer_results);
    priors_.push_back(one_layer_priors);
  }
  compFactors();
}

void CompProbValue::compFactors() {
  factors_.clear();
  // zero ptm
  int last_peak_index = peak_masses_.size() - 1;
  double hit_sum = 0;
  for (int i = 0; i < height_; i++) {
    hit_sum += results_[0][last_peak_index][i];
  }
  if (hit_sum == 0) {
    hit_sum = 1;
  }
  double zero_factor = 1.0 / hit_sum;
  factors_.push_back(zero_factor);

  double prior_sum = 0;
  for (int i = 0; i < height_; i++) {
    prior_sum += priors_[0][last_peak_index][i];
  }
  double peak_width = peak_mass_ends_[last_peak_index] 
      - peak_mass_bgns_[last_peak_index] + 1;
  for (int i = 1; i <= shift_num_; i++) {
    // i ptm
    double factor = 1.0 /prior_sum * i;
    factor = factor / pow(peak_width, i);
    factor = factor * K;
    factors_.push_back(factor);
  }
}

void CompProbValue::runClear(int page_pos) {
  std::memset(page_table_ + page_pos, 0, block_table_size_ * sizeof(double));
}

void CompProbValue::runFirstLayerInit(int win_table_bgn, int win_table_end) {
  double base_prob = 1.0;
  for (unsigned int i = 0; i < n_term_acid_masses_.size(); i++) {
    if (n_term_acid_frequencies_[i] <= 0) {
      continue;
    }
    //logger.debug("nTermAcidMass " + nTermAcidMasses[i] + " freq " +
    // nTermAcidFrequences[i]);
    int pos = n_term_acid_masses_[i];
    int k = pos * height_;
    // logger.debug("k " + k + " end " + blockTableSize);
    if (k >= win_table_bgn && k <= win_table_end) {
      page_table_[k % page_table_size_] += base_prob
          * n_term_acid_frequencies_[i];
      // logger.debug("pageTable value " + pageTable[k %
      // pageTableSize]);
    }
  }
}

void CompProbValue::runInit(std::vector<std::vector<double>> &results, 
                            int win_table_bgn, int win_table_end) {
  // zero
  double base_prob = 1.0;
  for (unsigned int i = 0; i < n_term_acid_masses_.size(); i++) {
    if (n_term_acid_frequencies_[i] <= 0) {
      continue;
    }
    int pos = n_term_acid_masses_[i];
    for (int k = pos * height_; k < pos * height_ + residue_avg_len_ * height_;
         k += height_) {
      if (k >= win_table_bgn && k <= win_table_end) {
        page_table_[k % page_table_size_] += base_prob * n_term_acid_frequencies_[i];
      }
      if (k > win_table_end) {
        break;
      }
    }
  }

  // results
  for (unsigned int i = 0; i < peak_masses_.size(); i++) {
    if (peak_table_ends_[i] >= win_table_end) {
      break;
    }
    if (peak_table_ends_[i] + residue_avg_len_ * height_ > win_table_bgn) {
      for (int k = peak_table_ends_[i] + 1; k < peak_table_ends_[i]
           + residue_avg_len_ * height_; k += height_) {
        if (k >= win_table_bgn && k <= win_table_end) {
          for (int h = 0; h < height_; h++) {
            page_table_[(k + h) % page_table_size_] += results[i][h];
          }
        }
        if (k > win_table_end) {
          break;
        }
      }
    }
  }
}

void CompProbValue::runAddProb(int page_pos, int prev_pos, int size, double f) {
  for (int i = page_pos; i < page_pos + size; i++) {
    page_table_[i] = page_table_[i] + page_table_[prev_pos] * f;
    prev_pos++;
  }
}

void CompProbValue::updateCol(int col_end, int scr) {
  int col_bgn = col_end - height_ + 1;
  int j, pre;
  for (j = 1; j <= scr; j++) {
    pre = col_end - j;
    if (pre >= col_bgn) {
      page_table_[col_end] += page_table_[pre];
    }
  }
  for (j = col_end - scr - 1; j >= col_bgn; j--) {
    page_table_[j + scr] = page_table_[j];
  }
  for (j = col_bgn; j < col_bgn + scr; j++) {
    if (j <= col_end) {
      page_table_[j] = 0;
    }
  }
}

void CompProbValue::runUpdate(int page_end, int scr_bgn, int scr_end) {
  int col_end = page_end;
  for (int i = scr_end; i >= scr_bgn; i--) {
    int scr = pos_scores_[i];
    if (scr > 0) {
      updateCol(col_end, scr);
    }
    col_end = col_end - height_;
  }
}

void CompProbValue::compOneLayer(std::vector<std::vector<double>> &prev_results, 
                                 bool is_first_layer,  
                                 std::vector<std::vector<double>> &cur_results, 
                                 std::vector<std::vector<double>> &cur_priors) {
  cur_results.clear();
  cur_priors.clear();
  for (unsigned int i = 0; i < peak_masses_.size(); i++) {
    std::vector<double> tmp(height_, 0.0);
    cur_results.push_back(tmp);
    cur_priors.push_back(tmp);
  }
  unsigned int peak_index = 0;
  for (int win_table_bgn = 0; win_table_bgn < sp_table_size_; 
       win_table_bgn = win_table_bgn + block_table_size_) {
    // clear result block
    int win_table_end = win_table_bgn + block_table_size_ - 1;
    int page_pos = win_table_bgn % page_table_size_;
    runClear(page_pos);
    if (is_first_layer) {
      runFirstLayerInit(win_table_bgn, win_table_end);
    } else {
      runInit(prev_results, win_table_bgn, win_table_end);
    }
    for (unsigned int i = 0; i < acid_dists_.size(); i++) {
      if (residue_frequencies_[i] <= 0) {
        continue;
      }
      int prev_block_bgn = page_pos - acid_dists_[i];
      if (prev_block_bgn < 0) {
        prev_block_bgn += page_table_size_;
      }
      int prev_block_end = prev_block_bgn + block_table_size_;
      if (prev_block_end < page_table_size_) {
        runAddProb(page_pos, prev_block_bgn, block_table_size_,
                   residue_frequencies_[i]);
      } else {
        int first_part_size = page_table_size_ - prev_block_bgn;
        runAddProb(page_pos, prev_block_bgn, first_part_size, 
                   residue_frequencies_[i]);
        runAddProb(page_pos + first_part_size, 0, block_table_size_
                   - first_part_size, residue_frequencies_[i]);
      }
    }
    int sp_bgn = win_table_bgn / height_;
    int sp_end = win_table_bgn / height_ + block_len_ - 1;
    runUpdate(page_pos + block_table_size_ - 1, sp_bgn, sp_end);
    // update peakPnt and get probs 
    while (peak_index < peak_masses_.size()
           && peak_table_ends_[peak_index] <= win_table_end) {
      for (int i = peak_mass_bgns_[peak_index]; i <= peak_mass_ends_[peak_index]; i++) {
        for (int j = 0; j < height_; j++) {
          int pos = (i * height_ + j) % page_table_size_;
          cur_results[peak_index][j] += page_table_[pos];
        }
      }

      int prior_mass_bgn = peak_mass_bgns_[peak_index] - residue_avg_len_;
      int prior_mass_end = peak_mass_bgns_[peak_index] - 1;
      // double sum = 0;
      for (int i = prior_mass_bgn; i <= prior_mass_end; i++) {
        for (int j = 0; j < height_; j++) {
          int pos = (i * height_ + j) % page_table_size_;
          if (pos < 0) {
            pos = pos + page_table_size_;
          }
          cur_priors[peak_index][j] += page_table_[pos];
        }
      }
      peak_index++;
    }
  }
}

double CompProbValue::getRawProb(int shift, int thresh) {
  double prob_sum = 0;
  if (results_.size() == 0) {
    return 1.0;
  }
  int last_peak_index = peak_masses_.size() - 1;
  for (int score = thresh; score < height_; score++) {
    prob_sum = prob_sum + results_[shift][last_peak_index][score];
  }
  return prob_sum * factors_[shift];
}

// The difference between getProb and getRawProb is nTermAcidFreqSum
double CompProbValue::getProb(int shift, int thresh) {
  if (results_.size() == 0) {
    return 1.0;
  }
  double raw_prob = getRawProb(shift, thresh);
  return raw_prob * n_term_acid_freq_sum_;
}

double CompProbValue::getOneValueProb(int shift, int value) {
  if (results_.size() == 0) {
    return 1.0;
  }
  int last_peak_index = peak_masses_.size() - 1;
  double prob = results_[shift][last_peak_index][value] * factors_[shift];
  return prob * n_term_acid_freq_sum_;
}

int getMaxScore(PrSMPtrVec prsms) {
  int score = 0;
  for (unsigned int i = 0; i < prsms.size(); i++) {
    if (prsms[i]->getMatchFragNum() > score) {
      score = (int) prsms[i]->getMatchFragNum();
    }
  }
  return score;
}

int getMaxShift(PrSMPtrVec prsms) {
  int shift = 0;
  for (unsigned int i = 0; i < prsms.size(); i++) {
    if (prsms[i]->getProteoformPtr()->getUnexpectedChangeNum() > shift) {
      shift = prsms[i]->getProteoformPtr()->getUnexpectedChangeNum();
    }
  }
  return shift;
}

void compProbArray(CompProbValue &comp_prob, PrmPeakPtrVec &peaks, 
                   PrSMPtrVec &prsms, bool strict, std::vector<double> &results) {
  int max_score = getMaxScore(prsms);
  int max_shift = getMaxShift(prsms);
  comp_prob.compute(peaks, max_score, max_shift, strict);
  results.clear();
  for (unsigned int i = 0; i < prsms.size(); i++) {
    int shift_num = prsms[i]->getProteoformPtr()->getUnexpectedChangeNum();
    int score = (int)prsms[i]->getMatchFragNum();
    results.push_back(comp_prob.getProb(shift_num, score));
  }
}

}
