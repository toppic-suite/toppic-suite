//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.


#include <cmath>

#include <cstring>

#include "base/logger.hpp"
#include "tdgf/comp_prob_value.hpp"

namespace prot {

ProbPeak::ProbPeak(PrmPeakPtr peak_ptr, int spectrum_id, int height, 
                   bool strict, double convert_ratio) {

  mass_ = (int) std::round(peak_ptr->getMonoMass() * convert_ratio);
  // we use NStrict and CRelax tolerance
  if (strict) {
    tolerance_ = std::ceil(peak_ptr->getStrictTolerance() * convert_ratio);
  } else {
    tolerance_ = std::ceil(peak_ptr->getNStrictCRelaxTolerance() * convert_ratio);
  }
  base_type_ptr_ = peak_ptr->getBaseTypePtr();
  spectrum_id_ = spectrum_id;
  mass_bgn_ = mass_ - tolerance_;
  mass_end_ = mass_ + tolerance_;
  table_bgn_ = mass_bgn_ * height;
  table_end_ = mass_end_ * height + (height - 1);
}

CompProbValue::CompProbValue(double convert_ratio, const ResFreqPtrVec &residue_ptrs, 
                             double residue_avg_len, int max_layer_num, int max_table_height, 
                             double max_sp_prec_mass) {
  convert_ratio_ = convert_ratio;
  residue_avg_len_ = residue_avg_len;

  for (size_t i = 0; i < residue_ptrs.size(); i++) {
    int int_mass = (int)std::round(residue_ptrs[i]->getMass() * convert_ratio);
    residue_masses_.push_back(int_mass);
    residue_frequencies_.push_back(residue_ptrs[i]->getFreq());
  }

  max_layer_num_ = max_layer_num;
  max_table_height_ = max_table_height;

  page_len_ = ORI_PAGE_LEN * (int) convert_ratio_;
  block_len_ = ORI_BLOCK_LEN * (int) convert_ratio_;
  page_table_ = new double[page_len_ * max_table_height_];

  max_sp_len_ = (int) std::round(max_sp_prec_mass * convert_ratio_);
  pos_scores_ = new short[max_sp_len_ + block_len_];
  tmp_pos_scores_ = new short[max_sp_len_ + block_len_];
  setFactors();
}

CompProbValue::~CompProbValue() {
  if (page_table_ != nullptr) {
    delete page_table_;
  }
  if (pos_scores_ != nullptr) {
    delete pos_scores_;
  }
  if (tmp_pos_scores_ != nullptr) {
    delete tmp_pos_scores_;
  }
}


inline void CompProbValue::setFactors () {
  //zero ptm 
  factors_.push_back(1);
  //one ptm
  factors_.push_back(K1());
  //i>=2 ptms based on experience: divide K by i^2
  int max_shift_num = max_layer_num_ - 1;
  for (int i = 2; i <= max_shift_num; i++) {
    factors_.push_back(K2());
  }
}

inline void CompProbValue::clearVar() {
  n_term_acid_masses_.clear();
  n_term_acid_frequencies_.clear();
  // we do not need to clear residue_masses and residue_frequencies

  prob_peaks_.clear();

  acid_dists_.clear();   // acidMass * height;

  results_.clear(); // nLayer peak number, height;

  // the prob that a randam protein has the same precursor to the spectrum 
  prec_probs_.clear(); 

  //page_table_ = new double[page_len_ * max_table_height_];
  memset(pos_scores_, 0, (max_sp_len_ + block_len_) * sizeof(short));
  memset(page_table_, 0, page_len_ * max_table_height_ * sizeof(double));
}

void CompProbValue::compute(const ResFreqPtrVec &n_residue_ptrs, 
                            const PrmPeakPtrVec2D &peak_ptr_2d, 
                            int thresh, int shift_num, bool strict) {
  //clear variables
  clearVar();
  for (size_t i = 0; i < n_residue_ptrs.size(); i++) {
    int int_mass = (int)std::round(n_residue_ptrs[i]->getMass() * convert_ratio_);
    n_term_acid_masses_.push_back(int_mass);
    n_term_acid_frequencies_.push_back(n_residue_ptrs[i]->getFreq());
  }
  // if score is less than 1, we do not compute p-value 
  if (thresh < 1) {
    return;
  }
  // initialize height 
  height_ = thresh + 1;
  if (height_ > max_table_height_) {
    height_ = max_table_height_;
  }

  setMassErr(peak_ptr_2d, strict);
  int group_spec_num = peak_ptr_2d.size();
  setPosScores(prob_peaks_, group_spec_num);
  setVars(thresh);
  if (sp_len_ + residue_avg_len_ >= max_sp_len_) {
    LOG_ERROR("Spectral precursor mass is too large ");
    return;
  }
  if (shift_num + 1 > max_layer_num_) {
    LOG_ERROR("Number of unexpected PTMs is too large " << shift_num);
    return;
  }
  shift_num_ = shift_num;
  comp();
}

inline void CompProbValue::setMassErr(const PrmPeakPtrVec2D &peak_ptr_2d, bool strict) {
  //LOG_DEBUG("MAX MASS " << peak_ptrs[peak_ptrs.size() - 1]->getMonoMass());
  // add mass 0;
  ProbPeak zero_peak (peak_ptr_2d[0][0], 0, height_, strict, convert_ratio_);
  prob_peaks_.push_back(zero_peak);

  for (size_t spec_id = 0; spec_id < peak_ptr_2d.size(); spec_id++) {
    // skip mass 0 and precursor mass 
    int peak_num = peak_ptr_2d[spec_id].size();
    for (int peak_id = 1; peak_id < peak_num-1; peak_id++) {
      ProbPeak cur_peak(peak_ptr_2d[spec_id][peak_id], spec_id, height_, strict, convert_ratio_);
      prob_peaks_.push_back(cur_peak);
    }
  }
  // precursor mass 
  int last_peak_idx = peak_ptr_2d[0].size() - 1;
  ProbPeak prec_peak(peak_ptr_2d[0][last_peak_idx], 0, height_, strict, convert_ratio_);
  prob_peaks_.push_back(prec_peak);
}

inline void CompProbValue::updatePosScores(const std::vector<ProbPeak> &prob_peaks, 
                                           int spectrum_id) {
  int len = max_sp_len_ + block_len_;
  memset(tmp_pos_scores_, 0, len * sizeof(short));
  //LOG_DEBUG("pos scr length " << len);
  //mass 0 and residue sum mass are not used for scoring 
  for (size_t i = 1; i < prob_peaks.size() -1; i++) {
    if (prob_peaks[i].spectrum_id_ == spectrum_id) {
      // here we use ceil/floor not round since each unit is very accurate
      int bgn = prob_peaks[i].mass_bgn_;
      if (bgn < 0) {
        bgn = 0;
      }
      int end = prob_peaks[i].mass_end_;
      if (end >= len) {
        end = len - 1;
      }
      //LOG_DEBUG("peak " << i << " bgn " << bgn << " end " << end);
      if (prob_peaks[i].base_type_ptr_ == BasePeakType::ORIGINAL) {
        // the peak is m
        for (int p = bgn; p <= end; p++) {
          if (tmp_pos_scores_[p] == 0) {
            tmp_pos_scores_[p] = 1;
          }
          if (tmp_pos_scores_[p] == 2) {
            tmp_pos_scores_[p] = 3;
          }
        }
      } else {
        for (int p = bgn; p <= end; p++) {
          if (tmp_pos_scores_[p] == 0) {
            tmp_pos_scores_[p] = 2;
          }
          if (tmp_pos_scores_[p] == 1) {
            tmp_pos_scores_[p] = 3;
          }
        }
      }
    }
  }

  for (int i = 0; i < len; i++) {
    if (tmp_pos_scores_[i] == 1 || tmp_pos_scores_[i] == 2) {
      pos_scores_[i] += 1;
    } else if (tmp_pos_scores_[i] == 3) {
      pos_scores_[i] += 2;
    }
  }
}

inline void CompProbValue::setPosScores(const std::vector<ProbPeak> &prob_peaks, 
                                        int group_spec_num) {
  for (int i = 0; i < group_spec_num; i++) {
    updatePosScores(prob_peaks, i);
  }
}

inline void CompProbValue::setVars(int thresh) {
  int last_peak_index = prob_peaks_.size() - 1;
  int max_peak_mass = prob_peaks_[last_peak_index].mass_end_;
  sp_len_ = max_peak_mass + residue_avg_len_ + 2;
  sp_table_size_ = sp_len_ * height_;
  page_table_size_ = page_len_ * height_;
  block_table_size_ = block_len_ * height_;
  //LOG_DEBUG("max mass " << max_peak_mass << " sp len " << sp_len_
  //          << " blk size " << block_table_size_ << " page size " << page_table_size_);
  //LOG_DEBUG("height " << height_);
  for (size_t i = 0; i < residue_masses_.size(); i++) {
    acid_dists_.push_back(residue_masses_[i] * height_);
  }
}

inline double CompProbValue::getShiftProb() { 
  std::vector<double> probs (height_, 0);
  size_t peak_index = prob_peaks_.size()-1;

  int prior_mass_bgn = prob_peaks_[peak_index].mass_bgn_ - residue_avg_len_;
  int prior_mass_end = prob_peaks_[peak_index].mass_bgn_ - 1;
  // double sum = 0;
  for (int i = prior_mass_bgn; i <= prior_mass_end; i++) {
    for (int j = 0; j < height_; j++) {
      int pos = (i * height_ + j) % page_table_size_;
      if (pos < 0) {
        pos = pos + page_table_size_;
      }

      probs[j] += page_table_[pos];
    }
  }
  double prob = 0;
  for (int i = 0; i < height_; i++) {
    prob += probs[i];
  }
  return prob;
}

void CompProbValue::comp() {
  std::vector<std::vector<double>> one_layer_results;
  //LOG_DEBUG("Start computation");
  std::vector<std::vector<double>> empty_vec;
  compOneLayer(empty_vec, true, one_layer_results);
  shift_prob_ = getShiftProb();
  //LOG_DEBUG("First level completed!");
  // score probability for each peak
  results_.push_back(one_layer_results);
  // prior probability for next layer
  for (int i = 1; i < shift_num_ + 1; i++) {
    one_layer_results.clear();
    compOneLayer(results_[i - 1], false, one_layer_results);
    results_.push_back(one_layer_results);
    //LOG_DEBUG("The " << i << " level completed!");
  }
  compPrecProbs();
  //LOG_DEBUG("compute precursor probabilities complete.");
}

void CompProbValue::compPrecProbs() {
  // zero ptm
  int last_peak_index = prob_peaks_.size() - 1;
  double prec_mass_minus_water = prob_peaks_[last_peak_index].mass_;
  double prob = 0;
  for (int i = 0; i < height_; i++) {
    prob += results_[0][last_peak_index][i];
  }
  if (prob == 0 && prec_mass_minus_water >= 1000) {
    LOG_WARN("Precursor probability is zero!"); 
  }
  prec_probs_.push_back(prob);

  // one ptm
  double peak_width = prob_peaks_[last_peak_index].mass_end_ 
      - prob_peaks_[last_peak_index].mass_bgn_ + 1;
  // because in results[], all probabilities are calculated 
  // with an error tolerance and all proteins are counted 
  // several times. In the prec_probability also needs to 
  // to be normalized by the error tolerance. 
  double one_ptm_prec_prob = shift_prob_ * peak_width;
  prec_probs_.push_back(one_ptm_prec_prob);

  // i >=2 ptms 
  for (int i = 2; i <= shift_num_; i++) {
    // for i ptms, all random proteins with an error 
    // < i * avg_residue_len are counted.
    double prec_prob = prob * i;
    // also, we need to normalize by peak_width ^ i;
    prec_prob = prec_prob * pow(peak_width, i);
    prec_probs_.push_back(prec_prob);
  }
}

inline void CompProbValue::runClear(int page_pos) {
  std::memset(&page_table_[page_pos], 0, block_table_size_ * sizeof(double));
}

inline void CompProbValue::runFirstLayerInit(int win_table_bgn, int win_table_end) {
  //LOG_DEBUG("N terminal acid masses number " << n_term_acid_masses_.size());
  for (size_t i = 0; i < n_term_acid_masses_.size(); i++) {
    if (n_term_acid_frequencies_[i] <= 0) {
      continue;
    }
    //LOG_DEBUG("nTermAcidMass " << n_term_acid_masses_[i]
    //          << " freq " << n_term_acid_frequencies_[i]);
    int pos = n_term_acid_masses_[i];
    int score = pos_scores_[pos];
    int k = pos * height_ + score;
    //LOG_DEBUG("k " << k << " end " << block_table_size_);
    if (k >= win_table_bgn && k <= win_table_end) {
      page_table_[k % page_table_size_] += n_term_acid_frequencies_[i];
      // logger.debug("pageTable value " + pageTable[k %
      // pageTableSize]);
    }
  }
}

inline void CompProbValue::runInit(std::vector<std::vector<double>> &results, 
                                   int win_table_bgn, int win_table_end) {
  // zero
  double base_prob = 1.0;
  for (size_t i = 0; i < n_term_acid_masses_.size(); i++) {
    if (n_term_acid_frequencies_[i] <= 0) {
      continue;
    }
    int pos = n_term_acid_masses_[i];
    for (int cur_pos = pos; cur_pos < pos + residue_avg_len_; cur_pos++) {
      int score = pos_scores_[cur_pos];
      int k = cur_pos * height_ + score;
      if (k >= win_table_bgn && k <= win_table_end) {
        page_table_[k % page_table_size_] += base_prob * n_term_acid_frequencies_[i];
      }
      if (k > win_table_end) {
        break;
      }
    }
  }

  // results
  for (size_t i = 0; i < prob_peaks_.size(); i++) {
    if (prob_peaks_[i].table_end_ >= win_table_end) {
      break;
    }
    if (prob_peaks_[i].table_end_ + residue_avg_len_ * height_ > win_table_bgn) {
      for (int k = prob_peaks_[i].table_end_ + 1; k < prob_peaks_[i].table_end_
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

inline void CompProbValue::runAddProb(int page_pos, int prev_pos, 
                                      int size, double f) {
  for (int i = page_pos; i < page_pos + size; i++) {
    page_table_[i] = page_table_[i] + page_table_[prev_pos] * f;
    prev_pos++;
  }
}

inline void CompProbValue::updateCol(int col_end, int scr) {
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

inline void CompProbValue::runUpdate(int page_end, int scr_bgn, int scr_end) {
  int col_end = page_end;
  for (int i = scr_end; i >= scr_bgn; i--) {
    int scr = pos_scores_[i];
    if (scr > 0) {
      updateCol(col_end, scr);
    }
    col_end = col_end - height_;
  }
}


inline void CompProbValue::compOneLayer(std::vector<std::vector<double>> &prev_results, 
                                        bool is_first_layer,  
                                        std::vector<std::vector<double>> &cur_results) { 
  cur_results.clear();
  //cur_priors.clear();
  for (size_t i = 0; i < prob_peaks_.size(); i++) {
    std::vector<double> tmp(height_, 0.0);
    cur_results.push_back(tmp);
  }
  //LOG_DEBUG("start for loop ");
  size_t peak_index = 0;
  for (int win_table_bgn = 0; win_table_bgn < sp_table_size_; 
       win_table_bgn = win_table_bgn + block_table_size_) {
    // clear result block
    int win_table_end = win_table_bgn + block_table_size_ - 1;
    int page_pos = win_table_bgn % page_table_size_;
    runClear(page_pos);
    //LOG_DEBUG("clear page complete ");
    if (is_first_layer) {
      runFirstLayerInit(win_table_bgn, win_table_end);
    } else {
      runInit(prev_results, win_table_bgn, win_table_end);
    }
    //LOG_DEBUG("unit page complete ");
    for (size_t i = 0; i < acid_dists_.size(); i++) {
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
    //LOG_DEBUG("add prob complete ");
    int sp_bgn = win_table_bgn / height_;
    int sp_end = win_table_bgn / height_ + block_len_ - 1;
    runUpdate(page_pos + block_table_size_ - 1, sp_bgn, sp_end);
    //LOG_DEBUG("update complete ");
    // update peakPnt and get probs 
    while (peak_index < prob_peaks_.size()
           && prob_peaks_[peak_index].table_end_ <= win_table_end) {
      for (int i = prob_peaks_[peak_index].mass_bgn_; i <= prob_peaks_[peak_index].mass_end_; i++) {
        for (int j = 0; j < height_; j++) {
          int pos = (i * height_ + j) % page_table_size_;
          cur_results[peak_index][j] += page_table_[pos];
        }
      }
      peak_index++;
    }
    //LOG_DEBUG("get prob complete ");
  }
}

double CompProbValue::getCondProb(int shift, int thresh) {
  double prob_sum = 0;
  if (results_.size() == 0) {
    return 1.0;
  }
  int last_peak_index = prob_peaks_.size() - 1;
  for (int score = thresh; score < height_; score++) {
    prob_sum = prob_sum + results_[shift][last_peak_index][score];
  }
  //compute conditional probability
  double cond_prob;
  if (prec_probs_[shift] > 0) {
    cond_prob = prob_sum / prec_probs_[shift];
  }
  else {
    cond_prob = 1000000;
  }
  //normalization
  double norm_cond_prob = cond_prob * factors_[shift];
  return norm_cond_prob;
}

double CompProbValue::getCondProbOneValue(int shift, int value) {
  if (results_.size() == 0) {
    return 1.0;
  }
  int last_peak_index = prob_peaks_.size() - 1;
  double prob = results_[shift][last_peak_index][value];
  //compute conditional probability
  double cond_prob;
  if (prec_probs_[shift] > 0) {
    cond_prob = prob / prec_probs_[shift];
  }
  else {
    cond_prob = 1000000;
  }
  //normalization
  double norm_cond_prob = cond_prob * factors_[shift];
  return norm_cond_prob;
}

int getMaxScore(const PrsmPtrVec &prsm_ptrs) {
  int score = 0;
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    if (prsm_ptrs[i]->getMatchFragNum() > score) {
      score = prsm_ptrs[i]->getMatchFragNum();
    }
  }
  return score;
}

int getMaxShift(PrsmPtrVec prsm_ptrs) {
  int shift = 0;
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    if (prsm_ptrs[i]->getProteoformPtr()->getMassShiftNum(MassShiftType::UNEXPECTED) > shift) {
      shift = prsm_ptrs[i]->getProteoformPtr()->getMassShiftNum(MassShiftType::UNEXPECTED);
    }
  }
  return shift;
}

void CompProbValue::compProbArray(CompProbValuePtr comp_prob_ptr, 
                                  const ResFreqPtrVec &n_term_residue_ptrs, 
                                  const PrmPeakPtrVec2D &peak_ptr_2d, 
                                  const PrsmPtrVec &prsm_ptrs, 
                                  bool strict, 
                                  std::vector<double> &results) {
  int max_score = getMaxScore(prsm_ptrs);
  int max_shift = getMaxShift(prsm_ptrs);
  //std::cout << std::endl << "max score " << max_score << " max shift " << max_shift << std::endl;
  comp_prob_ptr->compute(n_term_residue_ptrs, peak_ptr_2d, max_score, max_shift, strict);
  results.clear();
  for (size_t i = 0; i < prsm_ptrs.size(); i++) {
    int shift_num = prsm_ptrs[i]->getProteoformPtr()->getMassShiftNum(MassShiftType::UNEXPECTED);
    int score = prsm_ptrs[i]->getMatchFragNum();
    results.push_back(comp_prob_ptr->getCondProb(shift_num, score));
  }
}

}
