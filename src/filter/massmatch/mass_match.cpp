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

#include "filter/massmatch/mass_match.hpp"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>

#include "common/util/file_util.hpp"
#include "common/util/logger.hpp"

namespace toppic {


/* proteo_minus_water_masses[i]: the mass (minus water) of the ith proteoform
 * mass_2d[i]: a vector containing prefix residue masses of the ith proteoform
 * float_shift_2d[i]: a vector containing all possible shifts of the ith
 * proteoform pos_2d[i]: a vector containing the first residue position for each
 * shift.
 */
MassMatch::MassMatch(std::vector<double>& proteo_minus_water_masses,
                     std::vector<std::vector<int>>& mass_2d,
                     std::vector<std::vector<double>>& float_shift_2d,
                     std::vector<std::vector<int>>& pos_2d,
                     double max_proteoform_mass, double scale) {
  scale_ = scale;
  LOG_DEBUG("Scale: " << scale_);

  proteo_num_ = proteo_minus_water_masses.size();
  LOG_DEBUG("Proteoform number: " << proteo_minus_water_masses.size());
  proteo_minus_water_masses_ = proteo_minus_water_masses;

  col_num_ = max_proteoform_mass * scale_;
  LOG_DEBUG("column number: " << col_num_);

  LOG_DEBUG("start init");
  initProteoformBeginEnds(float_shift_2d);
  LOG_DEBUG("row number: " << row_num_);

  LOG_DEBUG("init indexes");
  initIndexes(mass_2d, float_shift_2d, pos_2d);
  initColIndexBounds();
}

void MassMatch::serializeMassMatch(std::string file_name,
                                   std::string dir_name) {
  std::string file_path = dir_name + file_util::getFileSeparator() + file_name;
  std::ofstream new_file(file_path, std::ofstream::binary);

  if (new_file.is_open()) {
    boost::archive::binary_oarchive oa(new_file, std::ios::binary);
    oa << scale_;
    oa << proteo_num_;
    oa << col_num_;
    oa << row_num_;

    oa << proteo_row_begins_;
    oa << proteo_row_ends_;
    oa << proteo_minus_water_masses_;
    oa << row_proteo_ids_;
    oa << trunc_shifts_;

    oa << col_index_begins_;
    oa << col_index_ends_;
    oa << col_indexes_;

    new_file.close();
  } else {
    LOG_ERROR("Failed to open the index file:" << file_path);
    exit(EXIT_FAILURE);
  }
}

void MassMatch::deserializeMassMatch(std::string new_file,
                                     std::string dir_name) {
  std::string file_path = dir_name + file_util::getFileSeparator() + new_file;
  std::ifstream file_to_read(file_path, std::ifstream::binary);

  if (file_to_read.is_open()) {
    boost::archive::binary_iarchive ia(file_to_read, std::ios::binary);

    ia >> scale_;
    ia >> proteo_num_;
    ia >> col_num_;
    ia >> row_num_;

    ia >> proteo_row_begins_;
    ia >> proteo_row_ends_;
    ia >> proteo_minus_water_masses_;
    ia >> row_proteo_ids_;
    ia >> trunc_shifts_;

    ia >> col_index_begins_;
    ia >> col_index_ends_;
    ia >> col_indexes_;

    file_to_read.close();
    initColIndexBounds();
  } else {
    LOG_ERROR("Failed to open the index file:" << file_path);
    exit(EXIT_FAILURE);
  }
}

void MassMatch::initProteoformBeginEnds(
    std::vector<std::vector<double>>& shift_2d) {
  // no need to init
  proteo_row_begins_.resize(proteo_num_);
  proteo_row_ends_.resize(proteo_num_);
  int pnt = 0;
  for (int i = 0; i < proteo_num_; i++) {
    proteo_row_begins_[i] = pnt;
    int len = shift_2d[i].size();
    proteo_row_ends_[i] = pnt + len - 1;
    pnt += len;
  }
  row_num_ = pnt;
  row_proteo_ids_.resize(row_num_);
  trunc_shifts_.resize(row_num_);
  for (int i = 0; i < proteo_num_; i++) {
    for (int j = proteo_row_begins_[i]; j <= proteo_row_ends_[i]; j++) {
      row_proteo_ids_[j] = i;
      int pos = j - proteo_row_begins_[i];
      trunc_shifts_[j] = shift_2d[i][pos];
    }
  }
}

inline std::vector<std::vector<int>> convertToInt(
    std::vector<std::vector<double>>& mass_2d, double scale) {
  std::vector<std::vector<int>> result;
  for (size_t i = 0; i < mass_2d.size(); i++) {
    std::vector<int> int_masses;
    for (size_t j = 0; j < mass_2d[i].size(); j++) {
      int value = std::floor(mass_2d[i][j] * scale + 0.5);
      int_masses.push_back(value);
    }
    result.push_back(int_masses);
  }
  return result;
}

void MassMatch::compColumnMatchNums(std::vector<std::vector<int>>& mass_2d,
                                    std::vector<std::vector<int>>& shift_2d,
                                    std::vector<std::vector<int>>& pos_2d,
                                    std::vector<int>& col_match_nums) {
  size_t proteo_num = mass_2d.size();
  for (size_t i = 0; i < proteo_num; i++) {
    for (size_t s = 0; s < shift_2d[i].size(); s++) {
      for (size_t cur = pos_2d[i][s]; cur < mass_2d[i].size(); cur++) {
        int shift_mass = mass_2d[i][cur] + shift_2d[i][s];
        if (shift_mass > 0) {
          if (shift_mass < col_num_) {
            col_match_nums[shift_mass]++;
          } else {
            break;
          }
        }
      }
    }
  }
}

void MassMatch::fillColumnIndex(std::vector<std::vector<int>>& mass_2d,
                                std::vector<std::vector<int>>& shift_2d,
                                std::vector<std::vector<int>>& pos_2d,
                                std::vector<int>& col_index_pnts) {
  for (size_t i = 0; i < mass_2d.size(); i++) {
    for (size_t s = 0; s < shift_2d[i].size(); s++) {
      for (size_t cur = pos_2d[i][s]; cur < mass_2d[i].size(); cur++) {
        int shift_mass = mass_2d[i][cur] + shift_2d[i][s];
        if (shift_mass > 0) {
          if (shift_mass < col_num_) {
            col_indexes_[col_index_pnts[shift_mass]] =
                proteo_row_begins_[i] + s;
            col_index_pnts[shift_mass]++;
          } else {
            break;
          }
        }
      }
    }
  }
}

void MassMatch::initIndexes(std::vector<std::vector<int>>& mass_2d,
                            std::vector<std::vector<double>>& float_shift_2d,
                            std::vector<std::vector<int>>& pos_2d) {
  std::vector<std::vector<int>> shift_2d = convertToInt(float_shift_2d, scale_);
  LOG_DEBUG("column num " << col_num_);
  std::vector<int> col_match_nums(col_num_, 0);
  // no need to initialize
  std::vector<int> col_index_pnts(col_num_);
  col_index_begins_.resize(col_num_);
  col_index_ends_.resize(col_num_);

  compColumnMatchNums(mass_2d, shift_2d, pos_2d, col_match_nums);

  int pnt = 0;
  for (int i = 0; i < col_num_; i++) {
    col_index_begins_[i] = pnt;
    col_index_pnts[i] = pnt;
    col_index_ends_[i] = pnt + col_match_nums[i] - 1;
    pnt += col_match_nums[i];
  }
  // no need to initialize
  col_indexes_.resize(pnt, 0);
  LOG_DEBUG("indexes size: " << pnt);
  fillColumnIndex(mass_2d, shift_2d, pos_2d, col_index_pnts);
}

void MassMatch::initColIndexBounds() {
  col_index_bounds_ = col_index_begins_;
  col_index_bounds_.push_back(static_cast<int>(col_indexes_.size()));
}

void MassMatch::compScores(
    const std::vector<std::pair<int, int>>& pref_mass_errors,
    std::vector<short>& scores) {
  compScores(pref_mass_errors, 0, 0.0, scores);
}

void MassMatch::compScores(
    const std::vector<std::pair<int, int>>& pref_mass_errors, int start,
    double shift, std::vector<short>& scores) {
  int begin_index;
  int end_index;
  int m;
  for (size_t i = start; i < pref_mass_errors.size(); i++) {
    m = pref_mass_errors[i].first + shift;
    // m - errors[i] performs better than m - errors[i] -  errors[bgn_pos]
    int left = m - pref_mass_errors[i].second;
    if (left < 0) {
      left = 0;
    }

    int right = m + pref_mass_errors[i].second;

    if (right < 0 || right >= col_num_) {
      continue;
    }

    begin_index = col_index_begins_[left];
    end_index = col_index_ends_[right];

    for (int j = begin_index; j <= end_index; j++) {
      scores[col_indexes_[j]]++;
    }
  }
}

// The anchor loop of DiagFilter::compute for the anchors [anchor_begin,
// anchor_end) with the score array of the current anchor in Score. The rows
// whose score reaches the threshold are noted by the increment loop itself
// (a rarely taken branch) and only they are folded into max_scores after the
// anchor, and the array is cleared with memset, so the per-anchor work on
// the whole array is one sequential write instead of a fill, a read and a
// read-modify-write of 16-bit scores. The column bounds and the first
// entries of the next pairs are prefetched: each pair reads them from memory
// the hardware prefetcher cannot predict.
template <typename Score>
void MassMatch::compDiagMaxScoresImpl(
    const std::vector<std::pair<int, int>>& pref_mass_errors, int anchor_begin,
    int anchor_end, int threshold, std::vector<Score>& scores,
    std::vector<int>& hot_rows, std::vector<short>& max_scores) {
  scores.assign(row_num_, 0);
  const int* bounds = col_index_bounds_.data();
  const int* entries = col_indexes_.data();
  Score* s = scores.data();
  short* mx = max_scores.data();
  const int n = static_cast<int>(pref_mass_errors.size());
  const int col_num = col_num_;
  const int row_num = row_num_;
  const int pf_bound = 8;  // pairs ahead for the column bound prefetch
  const int pf_entry = 3;  // pairs ahead for the entry prefetch
  const Score t = static_cast<Score>(threshold);
  for (int i = anchor_begin; i < anchor_end; i++) {
    if (i > anchor_begin) {
      std::memset(s, 0, sizeof(Score) * row_num);
    }
    hot_rows.clear();
    int shift = -pref_mass_errors[i].first;
    for (int k = i; k < n; k++) {
      if (k + pf_bound < n) {
        int m2 = pref_mass_errors[k + pf_bound].first + shift;
        int left2 = m2 - pref_mass_errors[k + pf_bound].second;
        int right2 = m2 + pref_mass_errors[k + pf_bound].second;
        if (left2 < 0) {
          left2 = 0;
        }
        if (right2 >= 0 && right2 < col_num) {
          __builtin_prefetch(bounds + left2);
          __builtin_prefetch(bounds + right2 + 1);
        }
      }
      if (k + pf_entry < n) {
        int m2 = pref_mass_errors[k + pf_entry].first + shift;
        int left2 = m2 - pref_mass_errors[k + pf_entry].second;
        int right2 = m2 + pref_mass_errors[k + pf_entry].second;
        if (left2 < 0) {
          left2 = 0;
        }
        if (right2 >= 0 && right2 < col_num) {
          int begin2 = bounds[left2];
          __builtin_prefetch(entries + begin2);
          __builtin_prefetch(entries + begin2 + 16);
        }
      }
      int m = pref_mass_errors[k].first + shift;
      int left = m - pref_mass_errors[k].second;
      if (left < 0) {
        left = 0;
      }
      int right = m + pref_mass_errors[k].second;
      if (right < 0 || right >= col_num) {
        continue;
      }
      int begin_index = bounds[left];
      int end_index = bounds[right + 1];
      for (int j = begin_index; j < end_index; j++) {
        int row = entries[j];
        if (++s[row] == t) {
          hot_rows.push_back(row);
        }
      }
    }
    for (size_t h = 0; h < hot_rows.size(); h++) {
      int row = hot_rows[h];
      if (s[row] > mx[row]) {
        mx[row] = static_cast<short>(s[row]);
      }
    }
  }
}

// An anchor's score for a row is at most the number of pairs of the anchor
// that hit the row, and a row is hit at most once by a pair when the column
// windows are narrower than the smallest residue mass (a row's entries are
// at least one residue mass apart). So the anchors with at most 255 later
// peaks are scored with 8-bit arrays and the earlier ones with 16-bit arrays.
void MassMatch::compDiagMaxScores(
    const std::vector<std::pair<int, int>>& pref_mass_errors, int threshold,
    std::vector<short>& max_scores) {
  // scratch arrays, kept per thread across spectrum sets and blocks
  static thread_local std::vector<uint8_t> scores8;
  static thread_local std::vector<uint16_t> scores16;
  static thread_local std::vector<int> hot_rows;
  max_scores.assign(row_num_, 0);
  const int n = static_cast<int>(pref_mass_errors.size());
  int max_error = 0;
  for (int k = 0; k < n; k++) {
    max_error = std::max(max_error, pref_mass_errors[k].second);
  }
  // glycine, 57.02 Da, scaled by the index scale
  const int min_residue_mass = static_cast<int>(57.0 * scale_);
  int first8 = n;
  if (2 * max_error + 1 < min_residue_mass) {
    first8 = std::max(0, n - 255);
  }
  if (first8 > 0) {
    compDiagMaxScoresImpl(pref_mass_errors, 0, first8, threshold, scores16,
                          hot_rows, max_scores);
  }
  if (first8 < n) {
    compDiagMaxScoresImpl(pref_mass_errors, first8, n, threshold, scores8,
                          hot_rows, max_scores);
  }
}

template <bool TrackPrecRows>
void MassMatch::updatePrecScoreImpl(const std::pair<int, int>& mass_error,
                                    std::vector<short>& scores,
                                    std::vector<int>* prec_rows) {
  // precursor mass
  int begin_index, end_index;
  int m = mass_error.first;
  // m - errors[i] performs better than m - errors[i] -  errors[bgn_pos]
  int left = m - mass_error.second;
  if (left < 0) {
    left = 0;
  }

  int right = m + mass_error.second;
  if (right >= col_num_) {
    right = col_num_ - 1;
  }

  if (left <= right) {
    // update scores
    begin_index = col_index_begins_[left];
    end_index = col_index_ends_[right];
    for (int j = begin_index; j <= end_index; j++) {
      int row = col_indexes_[j];
      if (scores[row] < getPrecursorMatchScore()) {
        scores[row] += getPrecursorMatchScore();
        if (TrackPrecRows) {
          prec_rows->push_back(row);
        }
      }
    }
  }
}

void MassMatch::updatePrecScore(const std::pair<int, int> mass_error,
                                std::vector<short>& scores) {
  updatePrecScoreImpl<false>(mass_error, scores, nullptr);
}

void MassMatch::compMatchScores(
    const std::vector<std::pair<int, int>>& pref_mass_errors,
    const std::pair<int, int>& prec_minus_water_mass_error,
    std::vector<short>& scores) {
  compScores(pref_mass_errors, 0, 0.0, scores);
  updatePrecScore(prec_minus_water_mass_error, scores);
}

void MassMatch::compMatchScores(
    const std::vector<std::pair<int, int>>& pref_mass_errors,
    const std::vector<std::pair<int, int>>& prec_minus_water_mass_errors,
    std::vector<short>& scores) {
  compScores(pref_mass_errors, 0, 0.0, scores);
  // precursor mass
  for (size_t i = 0; i < prec_minus_water_mass_errors.size(); i++) {
    std::pair<int, int> mass_error = prec_minus_water_mass_errors[i];
    updatePrecScore(mass_error, scores);
  }
}

void MassMatch::compMatchScores(
    const std::vector<std::pair<int, int>>& pref_mass_errors,
    const std::vector<std::pair<int, int>>& prec_minus_water_mass_errors,
    std::vector<short>& scores, std::vector<int>& prec_rows) {
  compScores(pref_mass_errors, 0, 0.0, scores);
  for (size_t i = 0; i < prec_minus_water_mass_errors.size(); i++) {
    updatePrecScoreImpl<true>(prec_minus_water_mass_errors[i], scores,
                              &prec_rows);
  }
}

} /* namespace toppic */
