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

#include "filter/massmatch/mass_match_util.hpp"

#include <algorithm>
#include <atomic>
#include <filesystem>
#include <functional>
#include <iostream>
#include <limits>

#include "common/util/file_util.hpp"
#include "common/util/logger.hpp"
#include "seq/prot_score.hpp"

namespace toppic {

namespace mass_match_util {

bool hasFilterLimits(const MassMatchPtr& index_ptr, int max_proteoform_mass,
                     int filter_scale) {
  return index_ptr->getScale() == filter_scale &&
         index_ptr->getColNum() == max_proteoform_mass * filter_scale;
}

void warnIndexLimits(const std::string& file_name,
                     const MassMatchPtr& index_ptr, int max_proteoform_mass,
                     int filter_scale) {
  static std::atomic<bool> warned{false};
  if (warned.exchange(true)) {
    return;
  }
  std::cout << "WARNING: The index file " << file_name << " has "
            << index_ptr->getColNum() << " columns at scale "
            << index_ptr->getScale() << ", not the "
            << max_proteoform_mass * filter_scale << " at scale "
            << filter_scale
            << " of this filter; the filter builds its indexes in memory. "
            << "Run topindex again." << std::endl;
}

void warnOldIndexLayout(const std::string& index_dir,
                        const std::string& suffix) {
  // the two names earlier versions gave the shared diagonal index
  static std::atomic<bool> warned{false};
  std::string dir = index_dir + file_util::getFileSeparator();
  if (!std::filesystem::exists(dir + "one_ptm_diag_index" + suffix) &&
      !std::filesystem::exists(dir + "multi_ptm_index" + suffix)) {
    return;
  }
  if (warned.exchange(true)) {
    return;
  }
  std::cout << "WARNING: The index directory " << index_dir
            << " was written by an earlier version (one_ptm_diag_index, "
            << "multi_ptm_index); the one-shift and multiple-shift filters "
            << "build their indexes in memory. Run topindex again, or rename "
            << "the one_ptm_diag_index files to ptm_diag_index and delete "
            << "the multi_ptm_index files." << std::endl;
}

// The best term_shift_num rows of a proteoform, ordered by score descending
// and, among equal scores, in the order in which the rows are visited. Only
// term_shift_num rows are kept, instead of sorting every row of the
// proteoform.
inline void addBestRow(int row, int score, int keep_num,
                       std::vector<std::pair<int, int>>& best) {
  if (static_cast<int>(best.size()) >= keep_num &&
      score <= best.back().second) {
    return;
  }
  size_t pos = best.size();
  while (pos > 0 && score > best[pos - 1].second) {
    pos--;
  }
  best.insert(best.begin() + pos, std::pair<int, int>(row, score));
  if (static_cast<int>(best.size()) > keep_num) {
    best.pop_back();
  }
}

void addTruncShifts(ProtCandidatePtrVec& prot_ptrs,
                    const MassMatchPtr& index_ptr, std::vector<short>& scores,
                    int term_shift_num, bool is_n_term) {
  const std::vector<int>& row_begins = index_ptr->getProteoRowBegins();
  const std::vector<int>& row_ends = index_ptr->getProteoRowEnds();
  const std::vector<double>& trunc_shifts = index_ptr->getTruncShifts();
  std::vector<std::pair<int, int>> best;
  for (size_t i = 0; i < prot_ptrs.size(); i++) {
    ProtCandidatePtr prot_ptr = prot_ptrs[i];
    int prot_id = prot_ptr->getProteinId();
    int bgn = row_begins[prot_id];
    int end = row_ends[prot_id];
    best.clear();
    for (int j = bgn; j <= end; j++) {
      addBestRow(j, scores[j], term_shift_num, best);
    }
    std::vector<double> shifts;
    for (size_t j = 0; j < best.size(); j++) {
      // check if the shift is the best score shift
      shifts.push_back(trunc_shifts[best[j].first]);
    }
    if (is_n_term) {
      prot_ptr->setNTermShifts(shifts);
    } else {
      prot_ptr->setCTermShifts(shifts);
    }
  }
}

class TruncScorePosition {
 public:
  TruncScorePosition(int score, int n_idx, int c_idx)
      : score_(score), n_idx_(n_idx), c_idx_(c_idx) {};

  int score_;
  int n_idx_;
  int c_idx_;
};

// The best term_group_num row pairs, ordered by score descending and, among
// equal scores, in the order in which the pairs are visited.
inline void addBestRowPair(const TruncScorePosition& pair_score, int keep_num,
                           std::vector<TruncScorePosition>& best) {
  if (static_cast<int>(best.size()) >= keep_num &&
      pair_score.score_ <= best.back().score_) {
    return;
  }
  size_t pos = best.size();
  while (pos > 0 && pair_score.score_ > best[pos - 1].score_) {
    pos--;
  }
  best.insert(best.begin() + pos, pair_score);
  if (static_cast<int>(best.size()) > keep_num) {
    best.pop_back();
  }
}

// The truncation shifts reported for the candidate proteoforms: the pairs of
// an N-terminal and a C-terminal truncation row whose proteoform mass is in
// the precursor mass window. Only the best term_group_num pairs are kept
// while the window is scanned, instead of collecting and sorting every pair.
std::vector<char> compNonincreasingShifts(const MassMatchPtr& index_ptr) {
  const std::vector<int>& row_begins = index_ptr->getProteoRowBegins();
  const std::vector<int>& row_ends = index_ptr->getProteoRowEnds();
  const std::vector<double>& trunc_shifts = index_ptr->getTruncShifts();
  std::vector<char> nonincr(row_begins.size(), 1);
  for (size_t i = 0; i < row_begins.size(); i++) {
    for (int j = row_begins[i] + 1; j <= row_ends[i]; j++) {
      if (trunc_shifts[j] > trunc_shifts[j - 1]) {
        nonincr[i] = 0;
        break;
      }
    }
  }
  return nonincr;
}

// The truncation window of a row, [k_end, k_start], is found by moving two
// pointers down from where the previous row left them: k_start to the first
// row of the reverse index whose mass reaches min_mass, k_end to the first
// whose mass exceeds max_mass. When the truncation shifts of the proteoform
// do not increase from row to row, the mass of a reverse row does not
// increase either (floating-point addition is monotone), so a pointer that
// passes over rows of the reverse index for a row would pass over them for
// every later row too: moving the pointers only for the rows whose window
// is scanned leaves them where moving them for every row does.
inline void moveTruncWindow(double proteo_mass, double n_trunc_shift,
                            const std::vector<double>& rev_trunc_shifts,
                            int rev_bgn, double min_mass, double max_mass,
                            int& k_start, int& k_end) {
  for (; k_start > rev_bgn; k_start--) {
    double mass = proteo_mass + n_trunc_shift + rev_trunc_shifts[k_start];
    if (mass >= min_mass) {
      break;
    }
  }
  for (; k_end > rev_bgn; k_end--) {
    double mass = proteo_mass + n_trunc_shift + rev_trunc_shifts[k_end];
    if (mass > max_mass) {
      break;
    }
  }
}

void addTruncShiftsWithRestriction(
    ProtCandidatePtrVec& prot_ptrs, const MassMatchPtr& index_ptr,
    const MassMatchPtr& rev_index_ptr,
    const std::vector<double>& proteo_minus_water_masses,
    std::vector<short>& scores, std::vector<short>& rev_scores,
    const std::vector<int>& rev_max_scores, double min_mass, double max_mass,
    int term_group_num, const std::vector<char>* nonincr_shifts) {
  const std::vector<int>& row_begins = index_ptr->getProteoRowBegins();
  const std::vector<int>& row_ends = index_ptr->getProteoRowEnds();
  const std::vector<double>& trunc_shifts = index_ptr->getTruncShifts();
  const std::vector<int>& rev_row_begins = rev_index_ptr->getProteoRowBegins();
  const std::vector<int>& rev_row_ends = rev_index_ptr->getProteoRowEnds();
  const std::vector<double>& rev_trunc_shifts = rev_index_ptr->getTruncShifts();

  int bgn, end, rev_bgn, rev_end;
  int k_start, k_end;
  std::vector<TruncScorePosition> best;
  for (size_t i = 0; i < prot_ptrs.size(); i++) {
    ProtCandidatePtr prot_ptr = prot_ptrs[i];
    int prot_id = prot_ptr->getProteinId();
    double proteo_mass = proteo_minus_water_masses[prot_id];
    bgn = row_begins[prot_id];
    end = row_ends[prot_id];
    rev_bgn = rev_row_begins[prot_id];
    rev_end = rev_row_ends[prot_id];

    k_start = rev_end;
    k_end = rev_end;

    best.clear();
    bool lazy = nonincr_shifts != nullptr && (*nonincr_shifts)[prot_id];
    for (int j = bgn; j <= end; j++) {
      // a row that cannot beat the worst kept pair adds nothing: once the
      // list is full, a row adds something when its score is at least
      // row_cutoff
      bool full = static_cast<int>(best.size()) >= term_group_num;
      int row_cutoff = full ? best.back().score_ + 1 - rev_max_scores[prot_id]
                            : std::numeric_limits<int>::min();
      if (lazy) {
        // the window is moved for the rows that are scanned only
        while (j <= end && scores[j] < row_cutoff) {
          j++;
        }
        if (j > end) {
          break;
        }
      }
      double n_trunc_shift = trunc_shifts[j];
      moveTruncWindow(proteo_mass, n_trunc_shift, rev_trunc_shifts, rev_bgn,
                      min_mass, max_mass, k_start, k_end);
      if (scores[j] < row_cutoff) {
        continue;
      }

      for (int k = k_start; k >= k_end; k--) {
        double c_trunc_shift = rev_trunc_shifts[k];
        double mass = proteo_mass + n_trunc_shift + c_trunc_shift;
        int score = scores[j] + rev_scores[k];
        if (mass >= min_mass && mass <= max_mass) {
          TruncScorePosition score_pos(score, j, k);
          addBestRowPair(score_pos, term_group_num, best);
        }
      }
    }
    // A row takes part in several of the best pairs - a C-terminal truncation
    // once for each N-terminal form, for example - and its shift is reported
    // once, in the order of its first pair (issue #436). The search drops a
    // repeated shift itself, so the result is the same.
    std::vector<double> n_shifts;
    std::vector<double> c_shifts;
    for (size_t j = 0; j < best.size(); j++) {
      bool new_n = true;
      bool new_c = true;
      for (size_t p = 0; p < j; p++) {
        new_n = new_n && best[p].n_idx_ != best[j].n_idx_;
        new_c = new_c && best[p].c_idx_ != best[j].c_idx_;
      }
      if (new_n) {
        n_shifts.push_back(trunc_shifts[best[j].n_idx_]);
      }
      if (new_c) {
        c_shifts.push_back(rev_trunc_shifts[best[j].c_idx_]);
      }
    }
    prot_ptr->setNTermShifts(n_shifts);
    prot_ptr->setCTermShifts(c_shifts);
  }
}

// The best keep_num proteoform scores, kept in a small ordered list while
// the proteoforms are scanned in id order: score descending and, among equal
// scores, protein id ascending - the order ProtCandidate::geneResults
// reports. The candidates are stored by value; a ProtScore object is only
// allocated for the survivors.
class BoundedProtScores {
 public:
  explicit BoundedProtScores(int keep_num) : keep_num_(keep_num) {}

  void add(int id, int score, double n_shift, double c_shift) {
    if (static_cast<int>(best_.size()) >= keep_num_ &&
        score <= best_.back().score_) {
      return;
    }
    size_t pos = best_.size();
    while (pos > 0 && score > best_[pos - 1].score_) {
      pos--;
    }
    best_.insert(best_.begin() + pos, Entry{id, score, n_shift, c_shift});
    if (static_cast<int>(best_.size()) > keep_num_) {
      best_.pop_back();
    }
  }

  ProtCandidatePtrVec toResults() const {
    ProtCandidatePtrVec results;
    for (size_t i = 0; i < best_.size(); i++) {
      ProtScorePtr prot_score = std::make_shared<ProtScore>(
          best_[i].id_, best_[i].score_, best_[i].n_shift_, best_[i].c_shift_);
      results.push_back(std::make_shared<ProtCandidate>(prot_score));
    }
    return results;
  }

 private:
  struct Entry {
    int id_;
    int score_;
    double n_shift_;
    double c_shift_;
  };
  int keep_num_;
  std::vector<Entry> best_;
};

ProtCandidatePtrVec findZeroShiftTopProteins(std::vector<short>& scores,
                                             std::vector<short>& rev_scores,
                                             const MassMatchPtr& index_ptr,
                                             const MassMatchPtr& rev_index_ptr,
                                             double prec_minus_water_mass,
                                             double prec_error_tole,
                                             int threshold, int top_num) {
  int prec_match_score = MassMatch::getPrecursorMatchScore();
  double min_mass = prec_minus_water_mass - prec_error_tole;
  double max_mass = prec_minus_water_mass + prec_error_tole;
  const std::vector<int>& row_begins = index_ptr->getProteoRowBegins();
  const std::vector<int>& row_ends = index_ptr->getProteoRowEnds();
  const std::vector<double>& proteo_minus_water_masses =
      index_ptr->getProteoMinusWaterMasses();
  const std::vector<double>& trunc_shifts = index_ptr->getTruncShifts();
  const std::vector<int>& rev_row_begins = rev_index_ptr->getProteoRowBegins();
  const std::vector<int>& rev_row_ends = rev_index_ptr->getProteoRowEnds();
  const std::vector<double>& rev_trunc_shifts = rev_index_ptr->getTruncShifts();

  ProtScorePtrVec proteo_scores;
  int bgn, end, rev_bgn, rev_end;
  std::vector<int> prec_match_rows;
  std::vector<int> rev_prec_match_rows;
  int best_score = 0;
  double best_n_trunc_shift = 0;
  double best_c_trunc_shift = 0;
  for (size_t i = 0; i < row_begins.size(); i++) {
    bgn = row_begins[i];
    end = row_ends[i];
    // get precursor match rows
    prec_match_rows.clear();
    for (int j = bgn; j <= end; j++) {
      if (scores[j] >= prec_match_score) {
        prec_match_rows.push_back(j);
      }
    }

    // get reverse precursor match rows
    rev_bgn = rev_row_begins[i];
    rev_end = rev_row_ends[i];
    rev_prec_match_rows.clear();
    for (int j = rev_bgn; j <= rev_end; j++) {
      if (rev_scores[j] >= prec_match_score) {
        rev_prec_match_rows.push_back(j);
      }
    }

    best_score = 0;
    best_n_trunc_shift = 0;
    best_c_trunc_shift = 0;
    for (size_t j = 0; j < prec_match_rows.size(); j++) {
      for (size_t k = 0; k < rev_prec_match_rows.size(); k++) {
        double n_trunc_shift = trunc_shifts[prec_match_rows[j]];
        double c_trunc_shift = rev_trunc_shifts[rev_prec_match_rows[k]];
        double mass =
            proteo_minus_water_masses[i] + n_trunc_shift + c_trunc_shift;
        double score =
            scores[prec_match_rows[j]] + rev_scores[rev_prec_match_rows[k]];
        if (mass >= min_mass && mass <= max_mass && score > best_score) {
          best_score = score;
          best_n_trunc_shift = n_trunc_shift;
          best_c_trunc_shift = c_trunc_shift;
        }
      }
    }
    if (best_score >= threshold) {
      ProtScorePtr proteo_score = std::make_shared<ProtScore>(
          i, best_score, best_n_trunc_shift, best_c_trunc_shift);
      proteo_scores.push_back(proteo_score);
    }
  }

  ProtCandidatePtrVec proteins =
      ProtCandidate::geneResults(proteo_scores, threshold, top_num);
  return proteins;
}

void compProteoMaxScores(const std::vector<short>& scores,
                         const MassMatchPtr& index_ptr,
                         std::vector<int>& max_scores) {
  const std::vector<int>& row_begins = index_ptr->getProteoRowBegins();
  const std::vector<int>& row_ends = index_ptr->getProteoRowEnds();
  max_scores.assign(row_begins.size(), 0);
  for (size_t i = 0; i < row_begins.size(); i++) {
    int best = 0;
    for (int j = row_begins[i]; j <= row_ends[i]; j++) {
      if (scores[j] > best) {
        best = scores[j];
      }
    }
    max_scores[i] = best;
  }
}

ProtCandidatePtrVec simpleFindZeroShiftTopProteins(
    const std::vector<int>& max_scores, const std::vector<int>& rev_max_scores,
    int threshold, int top_num) {
  BoundedProtScores best(top_num);
  // shift information is not used in simple filtering
  for (size_t i = 0; i < max_scores.size(); i++) {
    int best_total_score = max_scores[i] + rev_max_scores[i];
    if (best_total_score >= threshold) {
      best.add(i, best_total_score, 0.0, 0.0);
    }
  }
  return best.toResults();
}

inline ProtCandidatePtrVec findOneShiftTopProteinsWithoutRestrict(
    std::vector<short>& scores, std::vector<short>& rev_scores,
    const std::vector<int>& max_scores, const std::vector<int>& rev_max_scores,
    const MassMatchPtr& index_ptr, const MassMatchPtr& rev_index_ptr,
    int term_shift_num, int threshold, int top_num) {
  std::vector<std::pair<int, int>> proteo_scores;
  proteo_scores.reserve(max_scores.size());
  for (size_t i = 0; i < max_scores.size(); i++) {
    std::pair<int, int> proteo_score(i, max_scores[i] + rev_max_scores[i]);
    proteo_scores.push_back(proteo_score);
  }

  ProtCandidatePtrVec proteins =
      ProtCandidate::geneResults(proteo_scores, threshold, top_num);
  bool n_term = true;
  addTruncShifts(proteins, index_ptr, scores, term_shift_num, n_term);
  n_term = false;
  addTruncShifts(proteins, rev_index_ptr, rev_scores, term_shift_num, n_term);
  return proteins;
}

// A proteoform kept by the bounded candidate list of the one-shift search.
class ShiftCand {
 public:
  ShiftCand(int id, int score, double n_shift, double c_shift)
      : id_(id), score_(score), n_shift_(n_shift), c_shift_(c_shift) {}

  int id_;
  int score_;
  double n_shift_;
  double c_shift_;
};

// The order of ProtCandidate::geneResults: score descending, protein id
// ascending.
inline bool cmpShiftCand(const ShiftCand& a, const ShiftCand& b) {
  if (a.score_ != b.score_) {
    return a.score_ > b.score_;
  }
  return a.id_ < b.id_;
}

// The proteoform ids with a bound (best N-terminal row score + best
// C-terminal row score) of at least the threshold, ordered by decreasing
// bound. A counting sort is used because the bounds are small integers.
inline void sortProteoIdsByBound(const std::vector<int>& max_scores,
                                 const std::vector<int>& rev_max_scores,
                                 int threshold, std::vector<int>& order) {
  order.clear();
  int max_bound = threshold - 1;
  for (size_t i = 0; i < max_scores.size(); i++) {
    int bound = max_scores[i] + rev_max_scores[i];
    if (bound > max_bound) {
      max_bound = bound;
    }
  }
  if (max_bound < threshold) {
    return;
  }
  std::vector<int> bucket_begins(max_bound + 1, 0);
  for (size_t i = 0; i < max_scores.size(); i++) {
    int bound = max_scores[i] + rev_max_scores[i];
    if (bound >= threshold) {
      bucket_begins[bound]++;
    }
  }
  int pnt = 0;
  for (int bound = max_bound; bound >= threshold; bound--) {
    int cnt = bucket_begins[bound];
    bucket_begins[bound] = pnt;
    pnt += cnt;
  }
  order.resize(pnt);
  for (size_t i = 0; i < max_scores.size(); i++) {
    int bound = max_scores[i] + rev_max_scores[i];
    if (bound >= threshold) {
      order[bucket_begins[bound]] = static_cast<int>(i);
      bucket_begins[bound]++;
    }
  }
}

// The proteoforms are visited in decreasing order of their bound and the
// search stops when the bound can no longer reach the score of the last kept
// candidate, so that the truncation window scan below is run only for the
// proteoforms that can be reported. The bound is an upper bound of the total
// score of a proteoform, and the visiting order does not change the result
// because the kept candidates are ordered with the comparator of
// ProtCandidate::geneResults.
inline ProtCandidatePtrVec findOneShiftTopProteinsWithRestrict(
    std::vector<short>& scores, std::vector<short>& rev_scores,
    const std::vector<int>& max_scores, const std::vector<int>& rev_max_scores,
    const MassMatchPtr& index_ptr, const MassMatchPtr& rev_index_ptr,
    double prec_minus_water_mass, double prec_error_tole, double min_shift,
    double max_shift, int term_shift_num, int threshold, int top_num,
    const std::vector<char>* nonincr_shifts) {
  double min_mass = prec_minus_water_mass - max_shift - prec_error_tole;
  double max_mass = prec_minus_water_mass - min_shift + prec_error_tole;
  // LOG_DEBUG("min mass " << min_mass << " max mass " << max_mass);
  const std::vector<int>& row_begins = index_ptr->getProteoRowBegins();
  const std::vector<int>& row_ends = index_ptr->getProteoRowEnds();
  const std::vector<double>& proteo_minus_water_masses =
      index_ptr->getProteoMinusWaterMasses();
  const std::vector<double>& trunc_shifts = index_ptr->getTruncShifts();
  const std::vector<int>& rev_row_begins = rev_index_ptr->getProteoRowBegins();
  const std::vector<int>& rev_row_ends = rev_index_ptr->getProteoRowEnds();
  const std::vector<double>& rev_trunc_shifts = rev_index_ptr->getTruncShifts();

  std::vector<int> order;
  sortProteoIdsByBound(max_scores, rev_max_scores, threshold, order);

  std::vector<ShiftCand> cands;
  cands.reserve(top_num + 1);
  int bgn, end, rev_bgn, rev_end;
  int k_start, k_end;
  for (size_t idx = 0; idx < order.size(); idx++) {
    int i = order[idx];
    int rev_max_score = rev_max_scores[i];
    int bound = max_scores[i] + rev_max_score;
    int cutoff = static_cast<int>(cands.size()) >= top_num ? cands.back().score_
                                                           : threshold;
    if (bound < cutoff) {
      break;
    }

    double proteo_mass = proteo_minus_water_masses[i];
    bgn = row_begins[i];
    end = row_ends[i];
    rev_bgn = rev_row_begins[i];
    rev_end = rev_row_ends[i];
    k_start = rev_end;
    k_end = rev_end;
    int best_score = 0;
    double best_n_term_shift = 0;
    double best_c_term_shift = 0;
    bool lazy = nonincr_shifts != nullptr && (*nonincr_shifts)[i];
    for (int j = bgn; j <= end; j++) {
      // A row that cannot beat the best pair so far, or reach the cutoff,
      // adds nothing: a pair below the cutoff cannot be kept, and the pair
      // that gives the proteoform's best score is the first to reach it,
      // which no skipped row can be. A row adds something when its score
      // is at least row_cutoff.
      int row_cutoff = std::max(best_score + 1, cutoff) - rev_max_score;
      if (lazy) {
        // the window is moved for the rows that are scanned only, so the
        // rows before the next one to scan are passed over at once
        while (j <= end && scores[j] < row_cutoff) {
          j++;
        }
        if (j > end) {
          break;
        }
      }
      double n_trunc_shift = trunc_shifts[j];
      moveTruncWindow(proteo_mass, n_trunc_shift, rev_trunc_shifts, rev_bgn,
                      min_mass, max_mass, k_start, k_end);
      if (scores[j] < row_cutoff) {
        continue;
      }

      for (int k = k_start; k >= k_end; k--) {
        double c_trunc_shift = rev_trunc_shifts[k];
        double mass = proteo_mass + n_trunc_shift + c_trunc_shift;
        double score = scores[j] + rev_scores[k];
        if (mass >= min_mass && mass <= max_mass && score > best_score) {
          best_score = score;
          best_n_term_shift = n_trunc_shift;
          best_c_term_shift = c_trunc_shift;
        }
      }
    }

    if (best_score >= threshold) {
      ShiftCand cand(i, best_score, best_n_term_shift, best_c_term_shift);
      std::vector<ShiftCand>::iterator pos =
          std::upper_bound(cands.begin(), cands.end(), cand, cmpShiftCand);
      if (static_cast<int>(pos - cands.begin()) < top_num) {
        cands.insert(pos, cand);
        if (static_cast<int>(cands.size()) > top_num) {
          cands.pop_back();
        }
      }
    }
  }

  ProtCandidatePtrVec proteins;
  for (size_t i = 0; i < cands.size(); i++) {
    ProtScorePtr proteo_score = std::make_shared<ProtScore>(
        cands[i].id_, cands[i].score_, cands[i].n_shift_, cands[i].c_shift_);
    proteins.push_back(std::make_shared<ProtCandidate>(proteo_score));
  }
  // add trunc shifts
  addTruncShiftsWithRestriction(proteins, index_ptr, rev_index_ptr,
                                proteo_minus_water_masses, scores, rev_scores,
                                rev_max_scores, min_mass, max_mass,
                                term_shift_num, nonincr_shifts);
  return proteins;
}

ProtCandidatePtrVec findOneShiftTopProteins(
    std::vector<short>& scores, std::vector<short>& rev_scores,
    const std::vector<int>& max_scores, const std::vector<int>& rev_max_scores,
    const MassMatchPtr& index_ptr, const MassMatchPtr& rev_index_ptr,
    double prec_minus_water_mass, double prec_error_tole, double min_shift,
    double max_shift, int term_shift_num, int threshold, int top_num,
    const std::vector<char>* nonincr_shifts) {
  // if allowed mass shifts is small, use restricted search
  if ((max_shift - min_shift) <= 2000) {
    return findOneShiftTopProteinsWithRestrict(
        scores, rev_scores, max_scores, rev_max_scores, index_ptr,
        rev_index_ptr, prec_minus_water_mass, prec_error_tole, min_shift,
        max_shift, term_shift_num, threshold, top_num, nonincr_shifts);
  } else {
    ProtCandidatePtrVec cand_prots;
    double adjust_max_shift = 1000;
    if (adjust_max_shift > max_shift) {
      adjust_max_shift = max_shift;
    }
    double adjust_min_shift = -1000;
    if (adjust_min_shift < min_shift) {
      adjust_min_shift = min_shift;
    }
    if (adjust_max_shift > adjust_min_shift) {
      cand_prots = findOneShiftTopProteinsWithRestrict(
          scores, rev_scores, max_scores, rev_max_scores, index_ptr,
          rev_index_ptr, prec_minus_water_mass, prec_error_tole,
          adjust_min_shift, adjust_max_shift, term_shift_num, threshold,
          top_num, nonincr_shifts);
    }
    ProtCandidatePtrVec unrestrict_cand_prots =
        findOneShiftTopProteinsWithoutRestrict(
            scores, rev_scores, max_scores, rev_max_scores, index_ptr,
            rev_index_ptr, term_shift_num, threshold, top_num);
    // concatenate
    cand_prots.insert(cand_prots.end(), unrestrict_cand_prots.begin(),
                      unrestrict_cand_prots.end());
    return cand_prots;
  }
}

ProtCandidatePtrVec findVarPtmTopProteins(
    std::vector<short>& scores, std::vector<short>& rev_scores,
    const MassMatchPtr& index_ptr, const MassMatchPtr& rev_index_ptr,
    double prec_minus_water_mass, double prec_error_tole,
    std::vector<double>& ptm_shifts, int threshold, int top_num) {
  int prec_match_score = MassMatch::getPrecursorMatchScore();
  double min_mass = prec_minus_water_mass - prec_error_tole;
  double max_mass = prec_minus_water_mass + prec_error_tole;
  const std::vector<int>& row_begins = index_ptr->getProteoRowBegins();
  const std::vector<int>& row_ends = index_ptr->getProteoRowEnds();
  const std::vector<double>& proteo_minus_water_masses =
      index_ptr->getProteoMinusWaterMasses();
  const std::vector<double>& trunc_shifts = index_ptr->getTruncShifts();
  const std::vector<int>& rev_row_begins = rev_index_ptr->getProteoRowBegins();
  const std::vector<int>& rev_row_ends = rev_index_ptr->getProteoRowEnds();
  const std::vector<double>& rev_trunc_shifts = rev_index_ptr->getTruncShifts();

  ProtScorePtrVec proteo_scores;
  for (size_t i = 0; i < row_begins.size(); i++) {
    int bgn = row_begins[i];
    int end = row_ends[i];
    // get precursor match rows
    std::vector<int> prec_match_rows;
    for (int j = bgn; j <= end; j++) {
      if (scores[j] >= prec_match_score) {
        prec_match_rows.push_back(j);
      }
    }

    // get reverse precursor match rows
    int rev_bgn = rev_row_begins[i];
    int rev_end = rev_row_ends[i];
    std::vector<int> rev_prec_match_rows;
    for (int j = rev_bgn; j <= rev_end; j++) {
      if (rev_scores[j] >= prec_match_score) {
        rev_prec_match_rows.push_back(j);
      }
    }

    int best_score = 0;
    double best_n_trunc_shift = 0;
    double best_c_trunc_shift = 0;
    for (size_t j = 0; j < prec_match_rows.size(); j++) {
      for (size_t k = 0; k < rev_prec_match_rows.size(); k++) {
        double n_trunc_shift = trunc_shifts[prec_match_rows[j]];
        double c_trunc_shift = rev_trunc_shifts[rev_prec_match_rows[k]];
        double mass =
            proteo_minus_water_masses[i] + n_trunc_shift + c_trunc_shift;
        if (mass < 0) {
          continue;
        }
        double score =
            scores[prec_match_rows[j]] + rev_scores[rev_prec_match_rows[k]];
        bool mass_match = false;
        for (size_t s = 0; s < ptm_shifts.size(); s++) {
          double shift_mass = mass + ptm_shifts[s];
          // LOG_DEBUG("shift mass " << shift_mass << " min mass " << min_mass
          // << " max mass " << max_mass);
          if (shift_mass >= min_mass && shift_mass <= max_mass) {
            mass_match = true;
            break;
          }
        }
        if (mass_match && score > best_score) {
          best_score = score;
          best_n_trunc_shift = n_trunc_shift;
          best_c_trunc_shift = c_trunc_shift;
        }
      }
    }
    if (best_score >= threshold) {
      ProtScorePtr proteo_score = std::make_shared<ProtScore>(
          i, best_score, best_n_trunc_shift, best_c_trunc_shift);
      proteo_scores.push_back(proteo_score);
    }
  }

  ProtCandidatePtrVec proteins =
      ProtCandidate::geneResults(proteo_scores, threshold, top_num);
  return proteins;
}

// The noted rows are marked in a bit set of the index's rows and read back
// in ascending order, which also groups them by proteoform, and a count per
// proteoform gives the ranges. Putting each row into a list of its
// proteoform and sorting the lists cost 39 ns a row, a cache miss or two
// each, for 20,000 rows per spectrum and block with four PTMs (toppic_test
// report 0077_bfb05226). The row scan this replaced visited the rows of a
// proteoform in ascending order, and the search below keeps the first pair
// with the best score, so the order is part of the result.
void compPrecMatchRows(const MassMatchPtr& index_ptr,
                       const std::vector<int>& prec_rows,
                       PrecMatchRows& prec_match_rows) {
  const std::vector<int>& row_proteo_ids = index_ptr->getRowProteoIds();
  size_t proteo_num = index_ptr->getProteoRowBegins().size();
  size_t word_num = (row_proteo_ids.size() + 63) / 64;
  std::vector<uint64_t>& words = prec_match_rows.words_;
  words.assign(word_num, 0);
  for (size_t i = 0; i < prec_rows.size(); i++) {
    int row = prec_rows[i];
    words[row >> 6] |= uint64_t{1} << (row & 63);
  }
  std::vector<int>& rows = prec_match_rows.rows_;
  rows.clear();
  for (size_t w = 0; w < word_num; w++) {
    uint64_t word = words[w];
    while (word != 0) {
      rows.push_back(static_cast<int>(w * 64) + __builtin_ctzll(word));
      word &= word - 1;
    }
  }
  std::vector<int>& begins = prec_match_rows.begins_;
  begins.assign(proteo_num + 1, 0);
  for (size_t i = 0; i < rows.size(); i++) {
    begins[row_proteo_ids[rows[i]] + 1]++;
  }
  for (size_t i = 0; i < proteo_num; i++) {
    begins[i + 1] += begins[i];
  }
}

namespace {

// Whether mass plus one of the PTM shifts is inside [min_mass, max_mass]:
// the test of the pair loop below, and the one that decides every pair.
inline bool varPtmMassMatch(double mass, const std::vector<double>& ptm_shifts,
                            double min_mass, double max_mass) {
  for (size_t s = 0; s < ptm_shifts.size(); s++) {
    double shift_mass = mass + ptm_shifts[s];
    if (shift_mass >= min_mass && shift_mass <= max_mass) {
      return true;
    }
  }
  return false;
}

// The pairs of one proteoform are searched through the mass windows when the
// pair loop would test more than this many shift conditions.
const size_t kVarPtmWindowMinTests = 64;

// A mass window is widened by this much on both sides, so that it contains
// every mass the exact test accepts: the window ends are computed in a
// different order from that test, and their rounding errors are about 1e-9
// Da for the masses of the index.
const double kVarPtmWindowPad = 1e-4;

// A proteoform's rows of the second index, (truncation shift, position),
// in ascending order of shift. The shifts of a diagonal index's rows fall
// with the row (the n-th row of a proteoform removes its first n residues),
// so the reversed rows are usually sorted already.
inline void sortRevRows(const int* rev_rows, size_t rev_num,
                        const std::vector<double>& rev_trunc_shifts,
                        std::vector<std::pair<double, int>>& sorted_rev) {
  sorted_rev.clear();
  bool falling = true;
  for (size_t k = 1; k < rev_num && falling; k++) {
    falling = rev_trunc_shifts[rev_rows[k]] < rev_trunc_shifts[rev_rows[k - 1]];
  }
  if (falling) {
    for (size_t k = rev_num; k > 0; k--) {
      sorted_rev.emplace_back(rev_trunc_shifts[rev_rows[k - 1]],
                              static_cast<int>(k - 1));
    }
    return;
  }
  for (size_t k = 0; k < rev_num; k++) {
    sorted_rev.emplace_back(rev_trunc_shifts[rev_rows[k]],
                            static_cast<int>(k));
  }
  std::sort(sorted_rev.begin(), sorted_rev.end());
}

}  // namespace

// A pair of rows matches when the proteoform mass with the two truncation
// shifts plus one of the PTM shifts is within the precursor error tolerance,
// and each proteoform keeps the first pair, in the order of the rows, with
// the highest score. The pair loop tests every pair against every PTM shift,
// which with the 35 shifts of four PTMs was nearly all of the stage: almost
// no pair matches, so the loop runs through all the shifts for almost every
// pair (toppic_test report 0077_bfb05226). The masses that can match are the
// windows [min_mass - shift, max_mass - shift], one per PTM shift, which
// depend on the spectrum only and span a few hundred Da. For a row j of the
// first index the rows of the second index that put the mass inside that
// span are found by a binary search over their sorted shifts, and a row
// whose mass is inside one of the windows is decided by the loop's own test,
// in the order of the rows, with the loop's update rule, so the result is
// the same pair. A row j whose score cannot beat the best score so far is
// skipped, since only a strictly higher score replaces the best pair.
ProtCandidatePtrVec findVarPtmTopProteins(
    const std::vector<short>& scores, const std::vector<short>& rev_scores,
    const PrecMatchRows& prec_match_rows_2d,
    const PrecMatchRows& rev_prec_match_rows_2d,
    const MassMatchPtr& index_ptr, const MassMatchPtr& rev_index_ptr,
    double prec_minus_water_mass, double prec_error_tole,
    std::vector<double>& ptm_shifts, int threshold, int top_num) {
  double min_mass = prec_minus_water_mass - prec_error_tole;
  double max_mass = prec_minus_water_mass + prec_error_tole;
  const std::vector<double>& proteo_minus_water_masses =
      index_ptr->getProteoMinusWaterMasses();
  const std::vector<double>& trunc_shifts = index_ptr->getTruncShifts();
  const std::vector<double>& rev_trunc_shifts = rev_index_ptr->getTruncShifts();

  // the mass windows, padded, sorted and merged where they overlap
  std::vector<std::pair<double, double>> windows;
  for (size_t s = 0; s < ptm_shifts.size(); s++) {
    windows.emplace_back(min_mass - ptm_shifts[s] - kVarPtmWindowPad,
                         max_mass - ptm_shifts[s] + kVarPtmWindowPad);
  }
  std::sort(windows.begin(), windows.end());
  std::vector<double> win_lows;
  std::vector<double> win_highs;
  for (size_t w = 0; w < windows.size(); w++) {
    if (!win_highs.empty() && windows[w].first <= win_highs.back()) {
      win_highs.back() = std::max(win_highs.back(), windows[w].second);
    } else {
      win_lows.push_back(windows[w].first);
      win_highs.push_back(windows[w].second);
    }
  }
  // per proteoform: the second index's rows sorted by shift, and the
  // positions of the rows whose mass is inside a window
  std::vector<std::pair<double, int>> sorted_rev;
  std::vector<int> cand_pos;

  BoundedProtScores best(top_num);
  const std::vector<int>& begins = prec_match_rows_2d.begins_;
  const std::vector<int>& rev_begins = rev_prec_match_rows_2d.begins_;
  size_t proteo_num = begins.size() - 1;
  for (size_t i = 0; i < proteo_num; i++) {
    if (begins[i] == begins[i + 1] || rev_begins[i] == rev_begins[i + 1]) {
      continue;
    }
    // the rows of proteoform i in the two indexes, in ascending order
    const int* prec_match_rows = prec_match_rows_2d.rows_.data() + begins[i];
    size_t row_num = begins[i + 1] - begins[i];
    const int* rev_prec_match_rows =
        rev_prec_match_rows_2d.rows_.data() + rev_begins[i];
    size_t rev_num = rev_begins[i + 1] - rev_begins[i];
    int max_rev_score = 0;
    for (size_t k = 0; k < rev_num; k++) {
      int rev_score = rev_scores[rev_prec_match_rows[k]];
      max_rev_score = std::max(max_rev_score, rev_score);
    }
    bool use_windows = rev_num * ptm_shifts.size() > kVarPtmWindowMinTests;
    if (use_windows) {
      sortRevRows(rev_prec_match_rows, rev_num, rev_trunc_shifts, sorted_rev);
    }
    int best_score = 0;
    double best_n_trunc_shift = 0;
    double best_c_trunc_shift = 0;
    for (size_t j = 0; j < row_num; j++) {
      int score_j = scores[prec_match_rows[j]];
      if (score_j + max_rev_score <= best_score) {
        continue;
      }
      double n_trunc_shift = trunc_shifts[prec_match_rows[j]];
      double base_mass = proteo_minus_water_masses[i] + n_trunc_shift;
      if (!use_windows) {
        for (size_t k = 0; k < rev_num; k++) {
          int score = score_j + rev_scores[rev_prec_match_rows[k]];
          if (score <= best_score) {
            continue;
          }
          double c_trunc_shift = rev_trunc_shifts[rev_prec_match_rows[k]];
          double mass = base_mass + c_trunc_shift;
          if (mass < 0) {
            continue;
          }
          if (varPtmMassMatch(mass, ptm_shifts, min_mass, max_mass)) {
            best_score = score;
            best_n_trunc_shift = n_trunc_shift;
            best_c_trunc_shift = c_trunc_shift;
          }
        }
        continue;
      }
      // the rows whose mass is within the span of the windows, and of them
      // the ones inside a window
      double shift_low = win_lows.front() - base_mass;
      double shift_high = win_highs.back() - base_mass;
      auto it = std::lower_bound(sorted_rev.begin(), sorted_rev.end(),
                                 shift_low,
                                 [](const std::pair<double, int>& a, double v) {
                                   return a.first < v;
                                 });
      cand_pos.clear();
      for (; it != sorted_rev.end() && it->first <= shift_high; ++it) {
        double mass = base_mass + it->first;
        size_t w = std::upper_bound(win_lows.begin(), win_lows.end(), mass) -
                   win_lows.begin();
        if (w > 0 && mass <= win_highs[w - 1]) {
          cand_pos.push_back(it->second);
        }
      }
      if (cand_pos.empty()) {
        continue;
      }
      // the order of the rows
      std::sort(cand_pos.begin(), cand_pos.end());
      for (size_t c = 0; c < cand_pos.size(); c++) {
        int k = cand_pos[c];
        int score = score_j + rev_scores[rev_prec_match_rows[k]];
        if (score <= best_score) {
          continue;
        }
        double c_trunc_shift = rev_trunc_shifts[rev_prec_match_rows[k]];
        double mass = base_mass + c_trunc_shift;
        if (mass < 0) {
          continue;
        }
        if (varPtmMassMatch(mass, ptm_shifts, min_mass, max_mass)) {
          best_score = score;
          best_n_trunc_shift = n_trunc_shift;
          best_c_trunc_shift = c_trunc_shift;
        }
      }
    }
    if (best_score >= threshold) {
      best.add(i, best_score, best_n_trunc_shift, best_c_trunc_shift);
    }
  }

  return best.toResults();
}

// Find top proteins for diagonal filtering
ProtCandidatePtrVec findDiagTopProteins(
    const std::vector<short>& scores, const std::vector<int>& proteo_row_begins,
    const std::vector<int>& proteo_row_ends, int threshold, int top_num) {
  // The best top_num proteoforms are kept in a small ordered list while the
  // proteoforms are scanned, instead of collecting the score of every
  // proteoform and sorting them all. The kept order - score descending and,
  // among equal scores, protein id ascending, because the proteoforms are
  // visited in id order - is exactly the order ProtCandidate::geneResults
  // reports.
  std::vector<std::pair<int, int>> best;
  for (size_t i = 0; i < proteo_row_begins.size(); i++) {
    int bgn = proteo_row_begins[i];
    int end = proteo_row_ends[i];
    int best_score = 0;
    for (int j = bgn; j <= end; j++) {
      if (scores[j] > best_score) {
        best_score = scores[j];
      }
    }
    if (best_score >= threshold) {
      addBestRow(i, best_score, top_num, best);
    }
  }
  ProtCandidatePtrVec prot_results;
  for (size_t j = 0; j < best.size(); j++) {
    prot_results.push_back(
        std::make_shared<ProtCandidate>(best[j].first, best[j].second));
  }
  return prot_results;
}

}  // namespace mass_match_util

}  // namespace toppic
