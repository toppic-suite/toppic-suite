//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#include <algorithm>

#include "common/util/logger.hpp"
#include "seq/prot_score.hpp"
#include "filter/massmatch/mass_match_util.hpp"

namespace toppic {

namespace mass_match_util {

// Find top proteins for diagonal filtering
ProtCandidatePtrVec findTopProteins(std::vector<short> &scores,
                                    std::vector<int> &proteo_row_begins,
                                    std::vector<int> &proteo_row_ends,
                                    int threshold, int num) {
  std::vector<std::pair<int, int>> diag_scores;
  for (size_t i = 0; i < proteo_row_begins.size(); i++) {
    int bgn = proteo_row_begins[i];
    int end = proteo_row_ends[i];
    int best_score = 0;
    for (int j = bgn; j <= end; j++) {
      if (scores[j] > best_score) {
        best_score = scores[j];
      }
    }
    std::pair<int, int> diag_score(i, best_score);
    diag_scores.push_back(diag_score);
  }
  return ProtCandidate::geneResults(diag_scores, threshold, num);
}

inline bool cmpScore(const std::pair<int, int> &a, const std::pair<int, int> &b) {
  return a.second > b.second;
}

void addTruncShifts(ProtCandidatePtrVec &prot_ptrs, MassMatchPtr index_ptr,
                    std::vector<short> &scores, int term_shift_num, bool is_n_term) {
  std::vector<int> row_begins = index_ptr->getProteoRowBegins();
  std::vector<int> row_ends = index_ptr->getProteoRowEnds();
  std::vector<double> trunc_shifts = index_ptr->getTruncShifts();
  for (size_t i = 0; i < prot_ptrs.size(); i++) {
    ProtCandidatePtr prot_ptr = prot_ptrs[i];
    int prot_id = prot_ptr->getProteinId();
    int bgn = row_begins[prot_id];
    int end = row_ends[prot_id];
    std::vector<std::pair<int, int>> pos_scores;
    pos_scores.resize(end-bgn + 1);
    for (int j = bgn; j <= end; j++) {
      std::pair<int, int> cur_pos_score(j, scores[j]);
      pos_scores[j - bgn] = cur_pos_score;
    }
    std::sort(pos_scores.begin(), pos_scores.end(), cmpScore);
    std::vector<double> shifts;
    // get the best score shift;
    double best_score_shift; 
    if (is_n_term) {
      best_score_shift = prot_ptr->getNTermShifts()[0];
    }
    else {
      best_score_shift = prot_ptr->getCTermShifts()[0];
    }
    shifts.push_back(best_score_shift);
    for (size_t j = 0; j < pos_scores.size(); j++) {
      if (static_cast<int>(shifts.size()) >= term_shift_num) {
        break;
      }
      // check if the shift is the best score shift
      double shift = trunc_shifts[pos_scores[j].first];
      if (shift != best_score_shift) {
        shifts.push_back(shift); 
      }
    }
    if (is_n_term) {
      prot_ptr->setNTermShifts(shifts);
    } else {
      prot_ptr->setCTermShifts(shifts);
    }
  }
}

ProtCandidatePtrVec findTopProteins(std::vector<short> &scores,
                                    std::vector<short> &rev_scores,
                                    MassMatchPtr index_ptr,
                                    MassMatchPtr rev_index_ptr,
                                    double threshold, int num,
                                    bool add_shifts, int shift_num) {
  std::vector<int> row_begins = index_ptr->getProteoRowBegins();
  std::vector<int> row_ends = index_ptr->getProteoRowEnds();
  std::vector<int> rev_row_begins = rev_index_ptr->getProteoRowBegins();
  std::vector<int> rev_row_ends = rev_index_ptr->getProteoRowEnds();
  std::vector<std::pair<int, int> > proteo_scores;
  for (size_t i = 0; i < row_begins.size(); i++) {
    int bgn = row_begins[i];
    int end = row_ends[i];
    int best_score = 0;
    for (int j = bgn; j <= end; j++) {
      if (scores[j] > best_score) {
        best_score = scores[j];
      }
    }

    int rev_bgn = rev_row_begins[i];
    int rev_end = rev_row_ends[i];
    int rev_best_score = 0;
    for (int j = rev_bgn; j <= rev_end; j++) {
      if (rev_scores[j] > rev_best_score) {
        rev_best_score = rev_scores[j];
      }
    }

    std::pair<int, int> proteo_score(i, best_score + rev_best_score);
    proteo_scores.push_back(proteo_score);
  }

  ProtCandidatePtrVec proteins = ProtCandidate::geneResults(proteo_scores, threshold, num);
  if (add_shifts) {
    bool n_term = true;
    addTruncShifts(proteins, index_ptr, scores, shift_num, n_term);
    n_term = false;
    addTruncShifts(proteins, rev_index_ptr, rev_scores, shift_num, n_term);
  }
  return proteins;
}

ProtCandidatePtrVec findOneShiftTopProteins(std::vector<short> &scores,
                                            std::vector<short> &rev_scores,
                                            MassMatchPtr index_ptr,
                                            MassMatchPtr rev_index_ptr,
                                            double prec_minus_water_mass,
                                            double prec_error_tole,
                                            double min_shift, double max_shift,  
                                            int term_shift_num, 
                                            double threshold, int top_num) {

  double min_mass = prec_minus_water_mass - max_shift - prec_error_tole;
  double max_mass = prec_minus_water_mass - min_shift + prec_error_tole;
  std::vector<int> row_begins = index_ptr->getProteoRowBegins();
  std::vector<int> row_ends = index_ptr->getProteoRowEnds();
  std::vector<double> proteo_minus_water_masses = index_ptr->getProteoMinusWaterMasses();
  std::vector<double> trunc_shifts = index_ptr->getTruncShifts();
  std::vector<int> rev_row_begins = rev_index_ptr->getProteoRowBegins();
  std::vector<int> rev_row_ends = rev_index_ptr->getProteoRowEnds();
  std::vector<double> rev_trunc_shifts = rev_index_ptr->getTruncShifts();

  ProtScorePtrVec proteo_scores;
  for (size_t i = 0; i < row_begins.size(); i++) {
    int bgn = row_begins[i];
    int end = row_ends[i];
    int rev_bgn = rev_row_begins[i];
    int rev_end = rev_row_ends[i];

    int k_start = rev_bgn; 
    int k_end = rev_bgn;
    int best_score = 0;
    double best_n_term_shift = 0;
    double best_c_term_shift = 0;
    for (int j = bgn; j <= end; j++) {
      double n_trunc_shift = trunc_shifts[j];
      double proteo_mass = proteo_minus_water_masses[i];
      // find k_start 
      for (; k_start <= rev_end; k_start++) {
        double mass = proteo_mass + n_trunc_shift + rev_trunc_shifts[k_start]; 
        if (mass <= max_mass) {
          break;
        }
      }
      // find k_end
      for (; k_end <= rev_end; k_end++) {
        double mass = proteo_mass + n_trunc_shift + rev_trunc_shifts[k_end]; 
        if (mass < min_mass) {
          break;
        }
      }

      for (int k = k_start; k <= k_end; k++) {
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
      ProtScorePtr proteo_score = std::make_shared<ProtScore>(i, best_score, 
                                                              best_n_term_shift, 
                                                              best_c_term_shift);
      proteo_scores.push_back(proteo_score);
    }
  }

  ProtCandidatePtrVec proteins = ProtCandidate::geneResults(proteo_scores, threshold, top_num);
  // add trunc shifts
  bool is_n_term = true;
  addTruncShifts(proteins, index_ptr, scores, term_shift_num, is_n_term);
  is_n_term = false;
  addTruncShifts(proteins, rev_index_ptr, rev_scores, term_shift_num, is_n_term);
  return proteins;
}

ProtCandidatePtrVec findZeroShiftTopProteins(std::vector<short> &scores,
                                             std::vector<short> &rev_scores,
                                             MassMatchPtr index_ptr,
                                             MassMatchPtr rev_index_ptr,
                                             double prec_minus_water_mass,
                                             double prec_error_tole, 
                                             double threshold, int top_num) {

  int prec_match_score = MassMatch::getPrecursorMatchScore();
  double min_mass = prec_minus_water_mass - prec_error_tole;
  double max_mass = prec_minus_water_mass + prec_error_tole;
  std::vector<int> row_begins = index_ptr->getProteoRowBegins();
  std::vector<int> row_ends = index_ptr->getProteoRowEnds();
  std::vector<double> proteo_minus_water_masses = index_ptr->getProteoMinusWaterMasses();
  std::vector<double> trunc_shifts = index_ptr->getTruncShifts();
  std::vector<int> rev_row_begins = rev_index_ptr->getProteoRowBegins();
  std::vector<int> rev_row_ends = rev_index_ptr->getProteoRowEnds();
  std::vector<double> rev_trunc_shifts = rev_index_ptr->getTruncShifts();

  ProtScorePtrVec proteo_scores;
  for (size_t i = 0; i < row_begins.size(); i++) {
    int bgn = row_begins[i];
    int end = row_ends[i];
    // get precursor match rows
    std::vector<int> prec_match_rows;
    for (int j = bgn; j <= end; j++) {
      if (scores[j] > prec_match_score) {
        prec_match_rows.push_back(j);
      }
    }

    // get reverse precursor match rows
    int rev_bgn = rev_row_begins[i];
    int rev_end = rev_row_ends[i];
    std::vector<int> rev_prec_match_rows;
    for (int j = rev_bgn; j <= rev_end; j++) {
      if (rev_scores[j] > prec_match_score) {
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
        double mass = proteo_minus_water_masses[i] + n_trunc_shift + c_trunc_shift; 
        double score = scores[prec_match_rows[j]] + rev_scores[rev_prec_match_rows[k]];
        if (mass >= min_mass && mass <= max_mass && score > best_score) {
          best_score = score;
          best_n_trunc_shift = n_trunc_shift; 
          best_c_trunc_shift = c_trunc_shift; 
        }
      }
    }
    if (best_score >= threshold) {
      ProtScorePtr proteo_score = std::make_shared<ProtScore>(i, best_score, 
                                                              best_n_trunc_shift, 
                                                              best_c_trunc_shift);
      proteo_scores.push_back(proteo_score);
    }
  }

  ProtCandidatePtrVec proteins = ProtCandidate::geneResults(proteo_scores, threshold, top_num);
  return proteins;
}

} // namespace mass_match_util

} // namespace toppic 
