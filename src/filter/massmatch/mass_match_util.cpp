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
#include <iostream>

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
void addTruncShiftsWithRestriction(
    ProtCandidatePtrVec& prot_ptrs, const MassMatchPtr& index_ptr,
    const MassMatchPtr& rev_index_ptr,
    const std::vector<double>& proteo_minus_water_masses,
    std::vector<short>& scores, std::vector<short>& rev_scores,
    const std::vector<int>& rev_max_scores, double min_mass, double max_mass,
    int term_group_num) {
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
    for (int j = bgn; j <= end; j++) {
      double n_trunc_shift = trunc_shifts[j];
      // find k_start
      for (; k_start > rev_bgn; k_start--) {
        double mass = proteo_mass + n_trunc_shift + rev_trunc_shifts[k_start];
        if (mass >= min_mass) {
          break;
        }
      }
      // find k_end
      for (; k_end > rev_bgn; k_end--) {
        double mass = proteo_mass + n_trunc_shift + rev_trunc_shifts[k_end];
        if (mass > max_mass) {
          break;
        }
      }
      // a row that cannot beat the worst kept pair adds nothing
      if (static_cast<int>(best.size()) >= term_group_num &&
          scores[j] + rev_max_scores[prot_id] <= best.back().score_) {
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
    std::vector<double> n_shifts;
    std::vector<double> c_shifts;
    for (size_t j = 0; j < best.size(); j++) {
      n_shifts.push_back(trunc_shifts[best[j].n_idx_]);
      c_shifts.push_back(rev_trunc_shifts[best[j].c_idx_]);
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
    double max_shift, int term_shift_num, int threshold, int top_num) {
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
    for (int j = bgn; j <= end; j++) {
      double n_trunc_shift = trunc_shifts[j];
      // find k_start
      for (; k_start > rev_bgn; k_start--) {
        double mass = proteo_mass + n_trunc_shift + rev_trunc_shifts[k_start];
        if (mass >= min_mass) {
          break;
        }
      }
      // find k_end
      for (; k_end > rev_bgn; k_end--) {
        double mass = proteo_mass + n_trunc_shift + rev_trunc_shifts[k_end];
        if (mass > max_mass) {
          break;
        }
      }
      // A row that cannot beat the best pair so far, or reach the cutoff,
      // adds nothing: a pair below the cutoff cannot be kept, and the pair
      // that gives the proteoform's best score is the first to reach it,
      // which no skipped row can be.
      if (scores[j] + rev_max_score <= best_score ||
          scores[j] + rev_max_score < cutoff) {
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
                                term_shift_num);
  return proteins;
}

ProtCandidatePtrVec findOneShiftTopProteins(
    std::vector<short>& scores, std::vector<short>& rev_scores,
    const std::vector<int>& max_scores, const std::vector<int>& rev_max_scores,
    const MassMatchPtr& index_ptr, const MassMatchPtr& rev_index_ptr,
    double prec_minus_water_mass, double prec_error_tole, double min_shift,
    double max_shift, int term_shift_num, int threshold, int top_num) {
  // if allowed mass shifts is small, use restricted search
  if ((max_shift - min_shift) <= 2000) {
    return findOneShiftTopProteinsWithRestrict(
        scores, rev_scores, max_scores, rev_max_scores, index_ptr,
        rev_index_ptr, prec_minus_water_mass, prec_error_tole, min_shift,
        max_shift, term_shift_num, threshold, top_num);
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
          top_num);
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

void compPrecMatchRows(const MassMatchPtr& index_ptr,
                       const std::vector<int>& prec_rows,
                       std::vector<std::vector<int>>& prec_match_rows) {
  const std::vector<int>& row_proteo_ids = index_ptr->getRowProteoIds();
  size_t proteo_num = index_ptr->getProteoRowBegins().size();
  if (prec_match_rows.size() != proteo_num) {
    prec_match_rows.resize(proteo_num);
  }
  for (size_t i = 0; i < proteo_num; i++) {
    prec_match_rows[i].clear();
  }
  for (size_t i = 0; i < prec_rows.size(); i++) {
    prec_match_rows[row_proteo_ids[prec_rows[i]]].push_back(prec_rows[i]);
  }
  // the row scan this replaces visited the rows of a proteoform in ascending
  // order, and the search below keeps the first pair with the best score
  for (size_t i = 0; i < proteo_num; i++) {
    if (prec_match_rows[i].size() > 1) {
      std::sort(prec_match_rows[i].begin(), prec_match_rows[i].end());
    }
  }
}

ProtCandidatePtrVec findVarPtmTopProteins(
    const std::vector<short>& scores, const std::vector<short>& rev_scores,
    const std::vector<std::vector<int>>& prec_match_rows_2d,
    const std::vector<std::vector<int>>& rev_prec_match_rows_2d,
    const MassMatchPtr& index_ptr, const MassMatchPtr& rev_index_ptr,
    double prec_minus_water_mass, double prec_error_tole,
    std::vector<double>& ptm_shifts, int threshold, int top_num) {
  double min_mass = prec_minus_water_mass - prec_error_tole;
  double max_mass = prec_minus_water_mass + prec_error_tole;
  const std::vector<double>& proteo_minus_water_masses =
      index_ptr->getProteoMinusWaterMasses();
  const std::vector<double>& trunc_shifts = index_ptr->getTruncShifts();
  const std::vector<double>& rev_trunc_shifts = rev_index_ptr->getTruncShifts();

  BoundedProtScores best(top_num);
  for (size_t i = 0; i < prec_match_rows_2d.size(); i++) {
    const std::vector<int>& prec_match_rows = prec_match_rows_2d[i];
    const std::vector<int>& rev_prec_match_rows = rev_prec_match_rows_2d[i];
    if (prec_match_rows.empty() || rev_prec_match_rows.empty()) {
      continue;
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
