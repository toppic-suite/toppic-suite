// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "stat/tdgf/comp_prob_value.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

#include "common/util/logger.hpp"
#include "ms/factory/prm_peak_factory.hpp"

namespace toppic {

CompProbValue::CompProbValue(double convert_ratio,
                             const ResFreqPtrVec& residue_ptrs,
                             double residue_avg_len, int max_layer_num,
                             int max_table_height, double max_sp_prec_mass,
                             bool strict_complement) {
  convert_ratio_ = convert_ratio;
  strict_complement_ = strict_complement;
  setResidues(residue_ptrs, static_cast<int>(residue_avg_len));

  max_layer_num_ = max_layer_num;
  max_table_height_ = max_table_height;

  page_len_ = ORI_PAGE_LEN * static_cast<int>(convert_ratio_);
  block_len_ = ORI_BLOCK_LEN * static_cast<int>(convert_ratio_);

  max_sp_len_ = static_cast<int>(std::round(max_sp_prec_mass * convert_ratio_));
  pos_scores_ = new short[max_sp_len_ + block_len_];
  tmp_pos_scores_ = new short[max_sp_len_ + block_len_];
  setFactors();
}

void CompProbValue::setResidues(const ResFreqPtrVec& residue_ptrs,
                                int residue_avg_len) {
  residue_avg_len_ = residue_avg_len;
  residue_masses_.clear();
  residue_frequencies_.clear();
  for (size_t i = 0; i < residue_ptrs.size(); i++) {
    int int_mass = static_cast<int>(
        std::round(residue_ptrs[i]->getMass() * convert_ratio_));
    residue_masses_.push_back(int_mass);
    residue_frequencies_.push_back(residue_ptrs[i]->getFreq());
  }
}

CompProbValue::~CompProbValue() {
  if (pos_scores_ != nullptr) {
    delete[] pos_scores_;
  }
  if (tmp_pos_scores_ != nullptr) {
    delete[] tmp_pos_scores_;
  }
}

inline void CompProbValue::setFactors() {
  // zero ptm
  factors_.push_back(1);
  // one ptm
  factors_.push_back(K1());
  // i>=2 ptms based on experience: divide K by i^2
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

  results_.clear();  // nLayer peak number, height;

  // the prob that a randam protein has the same precursor to the spectrum
  prec_probs_.clear();

  // The score arrays and the page table are not cleared here: only the
  // region a computation reads is cleared, in compute(), once the peak
  // positions and the table height are known (clearUsedRegions).
}

void CompProbValue::compute(const ResFreqPtrVec& n_residue_ptrs,
                            const PrmPeakPtrVec2D& peak_ptr_2d, int thresh,
                            int shift_num, bool strict, double prec_mass,
                            const PeakTolerancePtr& tole_ptr) {
  // clear variables
  clearVar();
  // if peak number is 0, we do not compute p-value
  int peak_num = 0;
  for (size_t i = 0; i < peak_ptr_2d.size(); i++) {
    peak_num += peak_ptr_2d[i].size();
  }
  if (peak_num == 0) {
    return;
  }
  // if score is less than 1, we do not compute p-value
  if (thresh < 1) {
    return;
  }
  // init n-terminal residue frequencies
  for (size_t i = 0; i < n_residue_ptrs.size(); i++) {
    int int_mass =
        (int)std::round(n_residue_ptrs[i]->getMass() * convert_ratio_);
    n_term_acid_masses_.push_back(int_mass);
    n_term_acid_frequencies_.push_back(n_residue_ptrs[i]->getFreq());
  }
  // initialize height
  height_ = thresh + 1;
  if (height_ > max_table_height_) {
    height_ = max_table_height_;
  }

  setMassErr(peak_ptr_2d, strict, prec_mass, tole_ptr);
  // the extent of the position score arrays the computation touches: every
  // read is below sp_len_ (the last peak plus the average residue length)
  // rounded up to a block boundary
  int last_peak_index = prob_peaks_.size() - 1;
  pos_len_ = prob_peaks_[last_peak_index].mass_end_ + residue_avg_len_ +
             2 * block_len_;
  if (pos_len_ > max_sp_len_ + block_len_) {
    pos_len_ = max_sp_len_ + block_len_;
  }
  memset(pos_scores_, 0, pos_len_ * sizeof(short));
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
  setChunkVars();
  shift_num_ = shift_num;
  comp();
}

inline void CompProbValue::setMassErr(const PrmPeakPtrVec2D& peak_ptr_2d,
                                      bool strict, double prec_mass,
                                      const PeakTolerancePtr& tole_ptr) {
  PrmPeakPtr zero_peak_ptr =
      prm_peak_factory::getZeroPeakPtr(0, prec_mass, tole_ptr, 0.0);
  ProbPeak zero_peak(zero_peak_ptr, 0, height_, strict, strict_complement_,
                     convert_ratio_);
  prob_peaks_.push_back(zero_peak);

  for (size_t spec_id = 0; spec_id < peak_ptr_2d.size(); spec_id++) {
    for (size_t peak_id = 0; peak_id < peak_ptr_2d[spec_id].size(); peak_id++) {
      ProbPeak cur_peak(peak_ptr_2d[spec_id][peak_id], spec_id, height_, strict,
                        strict_complement_, convert_ratio_);
      prob_peaks_.push_back(cur_peak);
    }
  }
  // precursor mass
  PrmPeakPtr prec_peak_ptr =
      prm_peak_factory::getPrecPeakPtr(0, prec_mass, tole_ptr, 0.0);
  ProbPeak prec_peak(prec_peak_ptr, 0, height_, strict, strict_complement_,
                     convert_ratio_);
  prob_peaks_.push_back(prec_peak);
}

inline void CompProbValue::updatePosScores(
    const std::vector<ProbPeak>& prob_peaks, int spectrum_id) {
  int len = pos_len_;
  memset(tmp_pos_scores_, 0, len * sizeof(short));
  // mass 0 and residue sum mass are not used for scoring
  for (size_t i = 1; i < prob_peaks.size() - 1; i++) {
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

inline void CompProbValue::setPosScores(const std::vector<ProbPeak>& prob_peaks,
                                        int group_spec_num) {
  for (int i = 0; i < group_spec_num; i++) {
    updatePosScores(prob_peaks, i);
  }
}

inline void CompProbValue::setVars(int thresh) {
  int last_peak_index = prob_peaks_.size() - 1;
  int max_peak_mass = prob_peaks_[last_peak_index].mass_end_;
  sp_len_ = max_peak_mass + residue_avg_len_ + 2;
}

// The probability of the prior columns of the precursor peak (the
// residue_avg_len_ columns before it) in layer 0: per row the columns in
// order (summed during the layer-0 sweep), then the rows in order, as the
// whole-height sweep summed its page.
inline double CompProbValue::getShiftProb() {
  double prob = 0;
  for (int i = 0; i < height_; i++) {
    prob += shift_probs_[i];
  }
  return prob;
}

void CompProbValue::comp() {
  std::vector<std::vector<double>> one_layer_results;
  std::vector<std::vector<double>> empty_vec;
  compOneLayer(empty_vec, true, one_layer_results);
  shift_prob_ = getShiftProb();
  // score probability for each peak
  results_.push_back(one_layer_results);
  // prior probability for next layer
  for (int i = 1; i < shift_num_ + 1; i++) {
    one_layer_results.clear();
    compOneLayer(results_[i - 1], false, one_layer_results);
    results_.push_back(one_layer_results);
  }
  compPrecProbs();
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
    LOG_INFO("Precursor probability is zero!");
  }
  prec_probs_.push_back(prob);

  // one ptm
  double peak_width = prob_peaks_[last_peak_index].mass_end_ -
                      prob_peaks_[last_peak_index].mass_bgn_ + 1;
  // because in results[], all probabilities are calculated
  // with an error tolerance and all proteins are counted
  // several times. The prec_probability also needs to
  // to be normalized by the error tolerance.
  double one_ptm_prec_prob = shift_prob_ * peak_width;
  prec_probs_.push_back(one_ptm_prec_prob);

  // i >=2 ptms
  for (int i = 2; i <= shift_num_; i++) {
    // for i ptms, all random proteins with an error
    // < i * avg_residue_len are counted.
    double prec_prob = one_ptm_prec_prob * i;
    // because in results[], all probabilities are calculated
    // with an error tolerance and all proteins are counted
    // several times. The prec_probability also needs to
    // to be normalized by the error tolerance.
    prec_prob = prec_prob * peak_width;
    prec_probs_.push_back(prec_prob);
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
  // compute conditional probability
  double cond_prob;
  if (prec_probs_[shift] > 0) {
    cond_prob = prob_sum / prec_probs_[shift];
  } else {
    cond_prob = 1000000;
  }
  // normalization
  double norm_cond_prob = cond_prob * factors_[shift];
  return norm_cond_prob;
}

double CompProbValue::getCondProbOneValue(int shift, int value) {
  if (results_.size() == 0) {
    return 1.0;
  }
  int last_peak_index = prob_peaks_.size() - 1;
  double prob = results_[shift][last_peak_index][value];
  // compute conditional probability
  double cond_prob;
  if (prec_probs_[shift] > 0) {
    cond_prob = prob / prec_probs_[shift];
  } else {
    cond_prob = 1000000;
  }
  // normalization
  double norm_cond_prob = cond_prob * factors_[shift];
  return norm_cond_prob;
}

// The geometry of the sliced table (see the header)
void CompProbValue::setChunkVars() {
  int min_mass = block_len_;
  int max_mass = 1;
  for (size_t i = 0; i < residue_masses_.size(); i++) {
    if (residue_frequencies_[i] <= 0) {
      continue;
    }
    min_mass = std::min(min_mass, residue_masses_[i]);
    max_mass = std::max(max_mass, residue_masses_[i]);
  }
  chunk_block_cols_ = std::min(block_len_, min_mass);
  if (chunk_block_cols_ < 1) {
    chunk_block_cols_ = 1;
  }
  chunk_page_blocks_ = max_mass / chunk_block_cols_ + 3;
  // the whole-height sweep clears the tail of its circular page once per
  // computation, so a shift layer's first blocks read, through the wrap,
  // what the previous layer left in the last max_mass columns of the page:
  // kept here (the last pass of the previous layer over those page
  // columns) and put before column 0 of each shift layer's page
  chunk_tail_cols_ = max_mass;
  tail_capture_.assign(static_cast<size_t>(chunk_tail_cols_) * height_, 0.0);
  chunk_table_.assign(static_cast<size_t>(chunk_page_blocks_) *
                          chunk_block_cols_ * kChunkRows,
                      0.0);
  score_cols_.clear();
  max_pos_score_ = 0;
  for (int i = 0; i < pos_len_; i++) {
    if (pos_scores_[i] > 0) {
      score_cols_.push_back(i);
      max_pos_score_ =
          std::max(max_pos_score_, static_cast<int>(pos_scores_[i]));
    }
  }
  saved_rows_.assign(score_cols_.size() * std::max(max_pos_score_, 1), 0.0);
}

// One layer of the table: the slices of kChunkRows rows, each a sweep
// over the columns (sweepChunk). The per-peak results and, in layer 0, the
// prior-column sums of the precursor are accumulated slice by slice.
inline void CompProbValue::compOneLayer(
    std::vector<std::vector<double>>& prev_results, bool is_first_layer,
    std::vector<std::vector<double>>& cur_results) {
  cur_results.clear();
  for (size_t i = 0; i < prob_peaks_.size(); i++) {
    std::vector<double> tmp(height_, 0.0);
    cur_results.push_back(tmp);
  }
  if (is_first_layer) {
    shift_probs_.assign(height_, 0.0);
  }
  for (int s0 = 0; s0 < height_; s0 += kChunkRows) {
    int s1 = std::min(s0 + kChunkRows, height_);
    sweepChunk(prev_results, is_first_layer, cur_results, s0, s1);
  }
}

// The rows [s0, s1) of a layer's table over all columns, block by block in
// a circular page of chunk_page_blocks_ blocks of chunk_block_cols_
// columns, a column's rows contiguous (kChunkRows per column). For every
// block, in the order of the whole-height sweep: the block is zeroed, the
// initial probabilities are added (the N-terminal residues; in the shift
// layers also the previous layer's results spread over the columns after
// each peak), every residue adds its frequency times the cells one residue
// mass earlier (earlier blocks only), the columns with a position score
// move their values up by the score (the top row collecting what moves
// past it; the rows that come from the slice below were saved by that
// slice's sweep), and the peaks that end in the block collect their
// columns.
inline void CompProbValue::sweepChunk(
    std::vector<std::vector<double>>& prev_results, bool is_first_layer,
    std::vector<std::vector<double>>& cur_results, int s0, int s1) {
  const int rows = s1 - s0;
  const int block_cols = chunk_block_cols_;
  const int page_cols = chunk_page_blocks_ * block_cols;
  const int block_cells = block_cols * kChunkRows;
  const int page_cells = page_cols * kChunkRows;
  double* table = chunk_table_.data();
  std::fill(chunk_table_.begin(), chunk_table_.end(), 0.0);
  const int scr_max = std::max(max_pos_score_, 1);
  const int save_rows = s1 < height_ ? std::min(scr_max, rows) : 0;
  size_t peak_index = 0;
  size_t score_q = 0;
  const size_t prec_index = prob_peaks_.size() - 1;
  const int prior_bgn = prob_peaks_[prec_index].mass_bgn_ - residue_avg_len_;
  const int prior_end = prob_peaks_[prec_index].mass_bgn_ - 1;
  std::vector<double> col_new(kChunkRows);
  if (!is_first_layer) {
    for (int j = 0; j < chunk_tail_cols_; j++) {
      int c = j - chunk_tail_cols_;
      int slot = ((c % page_cols) + page_cols) % page_cols;
      double* dst = table + slot * kChunkRows;
      const double* src = &tail_capture_[static_cast<size_t>(j) * height_];
      for (int s = s0; s < s1; s++) {
        dst[s - s0] = src[s];
      }
    }
  }
  const int tail_page_bgn = page_len_ - chunk_tail_cols_;
  for (int c0 = 0; c0 < sp_len_; c0 += block_cols) {
    const int c1 = c0 + block_cols;
    const int block_idx = c0 / block_cols;
    double* block = table + (block_idx % chunk_page_blocks_) * block_cells;
    std::memset(block, 0, block_cells * sizeof(double));
    auto cell = [&](int c, int s) -> double& {
      return block[(c - c0) * kChunkRows + (s - s0)];
    };
    if (is_first_layer) {
      for (size_t i = 0; i < n_term_acid_masses_.size(); i++) {
        if (n_term_acid_frequencies_[i] <= 0) {
          continue;
        }
        int pos = n_term_acid_masses_[i];
        if (pos < c0 || pos >= c1) {
          continue;
        }
        int score = pos_scores_[pos];
        if (score >= s0 && score < s1) {
          cell(pos, score) += n_term_acid_frequencies_[i];
        }
      }
    } else {
      double base_prob = 1.0;
      for (size_t i = 0; i < n_term_acid_masses_.size(); i++) {
        if (n_term_acid_frequencies_[i] <= 0) {
          continue;
        }
        int pos = n_term_acid_masses_[i];
        int from = std::max(pos, c0);
        int to = std::min(pos + residue_avg_len_, c1);
        for (int cur_pos = from; cur_pos < to; cur_pos++) {
          int score = pos_scores_[cur_pos];
          if (score >= s0 && score < s1) {
            cell(cur_pos, score) += base_prob * n_term_acid_frequencies_[i];
          }
        }
      }
      for (size_t i = 0; i < prob_peaks_.size(); i++) {
        if (prob_peaks_[i].mass_end_ >= c1 - 1) {
          break;
        }
        int from = std::max(prob_peaks_[i].mass_end_ + 1, c0);
        int to = std::min(prob_peaks_[i].mass_end_ + residue_avg_len_ + 1, c1);
        for (int c = from; c < to; c++) {
          for (int h = s0; h < s1; h++) {
            cell(c, h) += prev_results[i][h];
          }
        }
      }
    }
    for (size_t r = 0; r < residue_masses_.size(); r++) {
      if (residue_frequencies_[r] <= 0) {
        continue;
      }
      const double f = residue_frequencies_[r];
      int src_col = c0 - residue_masses_[r];
      int src_slot = ((src_col % page_cols) + page_cols) % page_cols;
      int src = src_slot * kChunkRows;
      int first = std::min(block_cells, page_cells - src);
      for (int k = 0; k < first; k++) {
        block[k] = block[k] + table[src + k] * f;
      }
      for (int k = first; k < block_cells; k++) {
        block[k] = block[k] + table[k - first] * f;
      }
    }
    while (score_q < score_cols_.size() && score_cols_[score_q] < c0) {
      score_q++;
    }
    for (; score_q < score_cols_.size() && score_cols_[score_q] < c1;
         score_q++) {
      const int c = score_cols_[score_q];
      const int scr = pos_scores_[c];
      double* col = &cell(c, s0);
      const double* saved = &saved_rows_[score_q * scr_max];
      auto value = [&](int s) -> double {
        return s >= s0 ? col[s - s0] : saved[s - (s0 - scr_max)];
      };
      for (int s = s0; s < s1; s++) {
        double v;
        if (s == height_ - 1) {
          v = col[s - s0];
          for (int j = 1; j <= scr; j++) {
            int pre = s - j;
            if (pre >= 0) {
              v += value(pre);
            }
          }
        } else if (s < scr) {
          v = 0;
        } else {
          v = value(s - scr);
        }
        col_new[s - s0] = v;
      }
      if (save_rows > 0) {
        double* save = &saved_rows_[score_q * scr_max];
        for (int t = 0; t < scr_max; t++) {
          int s = s1 - scr_max + t;
          save[t] = s >= s0 ? col[s - s0] : 0.0;
        }
      }
      for (int s = s0; s < s1; s++) {
        col[s - s0] = col_new[s - s0];
      }
    }
    while (peak_index < prob_peaks_.size() &&
           prob_peaks_[peak_index].mass_end_ <= c1 - 1) {
      for (int i = prob_peaks_[peak_index].mass_bgn_;
           i <= prob_peaks_[peak_index].mass_end_; i++) {
        int slot = ((i % page_cols) + page_cols) % page_cols;
        const double* src = table + slot * kChunkRows;
        for (int s = s0; s < s1; s++) {
          cur_results[peak_index][s] += src[s - s0];
        }
      }
      peak_index++;
    }
    // the columns of this block that lie in the tail of the whole-height
    // sweep's page (every pass; the last one stays), for the next layer
    {
      int p0 = c0 % page_len_;
      for (int c = c0; c < c1; c++) {
        int p = p0 + (c - c0);
        if (p >= page_len_) {
          p -= page_len_;
        }
        if (p >= tail_page_bgn) {
          double* dst = &tail_capture_[static_cast<size_t>(p - tail_page_bgn) * height_];
          const double* src = &cell(c, s0);
          for (int s = s0; s < s1; s++) {
            dst[s] = src[s - s0];
          }
        }
      }
    }
    if (is_first_layer) {
      int from = std::max(prior_bgn, c0);
      int to = std::min(prior_end, c1 - 1);
      for (int i = from; i <= to; i++) {
        int slot = ((i % page_cols) + page_cols) % page_cols;
        const double* src = table + slot * kChunkRows;
        for (int s = s0; s < s1; s++) {
          shift_probs_[s] += src[s - s0];
        }
      }
    }
  }
}

}  // namespace toppic
