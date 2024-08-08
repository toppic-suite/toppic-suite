//Copyright (c) 2014 - 2023, The Trustees of Indiana University.
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
    for (size_t j = 0; j < pos_scores.size(); j++) {
      if (static_cast<int>(shifts.size()) >= term_shift_num) {
        break;
      }
      // check if the shift is the best score shift
      double shift = trunc_shifts[pos_scores[j].first];
      shifts.push_back(shift); 
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
  TruncScorePosition(double score, int n_idx, int c_idx):
    score_(score),
    n_idx_(n_idx),
    c_idx_(c_idx) {};
                     
  double score_;
  int n_idx_;
  int c_idx_;
};

inline bool cmpScoreIndex(const TruncScorePosition &a, const TruncScorePosition &b) {
  return a.score_ > b.score_;
}


void addTruncShiftsWithRestriction(ProtCandidatePtrVec &prot_ptrs, 
                                   MassMatchPtr index_ptr, MassMatchPtr rev_index_ptr,
                                   std::vector<double> &proteo_minus_water_masses,
                                   std::vector<short> &scores, std::vector<short> &rev_scores,
                                   double min_mass, double max_mass, 
                                   int term_group_num) {
  std::vector<int> row_begins = index_ptr->getProteoRowBegins();
  std::vector<int> row_ends = index_ptr->getProteoRowEnds();
  std::vector<double> trunc_shifts = index_ptr->getTruncShifts();
  std::vector<int> rev_row_begins = rev_index_ptr->getProteoRowBegins();
  std::vector<int> rev_row_ends = rev_index_ptr->getProteoRowEnds();
  std::vector<double> rev_trunc_shifts = rev_index_ptr->getTruncShifts();

  int bgn, end, rev_bgn, rev_end;
  int k_start, k_end; 
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
  
    std::vector<TruncScorePosition> score_idxes;
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

      for (int k = k_start; k >= k_end; k--) {
        double c_trunc_shift = rev_trunc_shifts[k];
        double mass = proteo_mass + n_trunc_shift + c_trunc_shift; 
        double score = scores[j] + rev_scores[k];
        if (mass >= min_mass && mass <= max_mass) {
          TruncScorePosition score_pos(score, j, k);
          score_idxes.push_back(score_pos);
        }
      }
    }
    std::sort(score_idxes.begin(), score_idxes.end(), cmpScoreIndex);
    std::vector<double> n_shifts;
    std::vector<double> c_shifts;
    for (size_t j = 0; j < score_idxes.size(); j++) {
      if (static_cast<int>(n_shifts.size()) >= term_group_num) {
        break;
      }
      double n_shift = trunc_shifts[score_idxes[j].n_idx_];
      n_shifts.push_back(n_shift); 
      double c_shift = rev_trunc_shifts[score_idxes[j].c_idx_];
      c_shifts.push_back(c_shift);
    }
    prot_ptr->setNTermShifts(n_shifts);
    prot_ptr->setCTermShifts(c_shifts);
  }
}
    
ProtCandidatePtrVec findZeroShiftTopProteins(std::vector<short> &scores,
                                             std::vector<short> &rev_scores,
                                             MassMatchPtr index_ptr,
                                             MassMatchPtr rev_index_ptr,
                                             double prec_minus_water_mass,
                                             double prec_error_tole, 
                                             int threshold, int top_num) {

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

ProtCandidatePtrVec simpleFindZeroShiftTopProteins(std::vector<short> &scores,
                                                   std::vector<short> &rev_scores,
                                                   MassMatchPtr index_ptr,
                                                   MassMatchPtr rev_index_ptr,
                                                   double prec_minus_water_mass,
                                                   double prec_error_tole, 
                                                   int threshold, int top_num) {

  std::vector<int> row_begins = index_ptr->getProteoRowBegins();
  std::vector<int> row_ends = index_ptr->getProteoRowEnds();
  std::vector<int> rev_row_begins = rev_index_ptr->getProteoRowBegins();
  std::vector<int> rev_row_ends = rev_index_ptr->getProteoRowEnds();

  ProtScorePtrVec proteo_scores;
  int bgn, end, rev_bgn, rev_end;
  int best_n_score, best_c_score, best_total_score;
  // shift information is not used in simple filtering
  double best_n_trunc_shift = 0.0;
  double best_c_trunc_shift = 0.0;
  for (size_t i = 0; i < row_begins.size(); i++) {
    bgn = row_begins[i];
    end = row_ends[i];
    // get precursor match rows
    best_n_score = 0;
    for (int j = bgn; j <= end; j++) {
      if (scores[j] > best_n_score)  {
        best_n_score = scores[j];
      }
    }

    // get reverse precursor match rows
    rev_bgn = rev_row_begins[i];
    rev_end = rev_row_ends[i];
    best_c_score = 0;
    for (int j = rev_bgn; j <= rev_end; j++) {
      if (rev_scores[j] > best_c_score) {
        best_c_score = rev_scores[j];
      }
    }

    best_total_score = best_n_score + best_c_score;
    if (best_total_score >= threshold) {
      ProtScorePtr proteo_score = std::make_shared<ProtScore>(i, best_total_score, 
                                                              best_n_trunc_shift, 
                                                              best_c_trunc_shift);
      proteo_scores.push_back(proteo_score);
    }
  }

  ProtCandidatePtrVec proteins = ProtCandidate::geneResults(proteo_scores, threshold, top_num);
  return proteins;
}

inline ProtCandidatePtrVec findOneShiftTopProteinsWithoutRestrict(std::vector<short> &scores,
                                                                  std::vector<short> &rev_scores,
                                                                  MassMatchPtr index_ptr,
                                                                  MassMatchPtr rev_index_ptr,
                                                                  int term_shift_num, 
                                                                  int threshold, int top_num) {
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

  ProtCandidatePtrVec proteins = ProtCandidate::geneResults(proteo_scores, threshold, top_num);
  bool n_term = true;
  addTruncShifts(proteins, index_ptr, scores, term_shift_num, n_term);
  n_term = false;
  addTruncShifts(proteins, rev_index_ptr, rev_scores, term_shift_num, n_term);
  return proteins;
}


inline ProtCandidatePtrVec findOneShiftTopProteinsWithRestrict(std::vector<short> &scores,
                                                               std::vector<short> &rev_scores,
                                                               MassMatchPtr index_ptr,
                                                               MassMatchPtr rev_index_ptr,
                                                               double prec_minus_water_mass,
                                                               double prec_error_tole,
                                                               double min_shift, double max_shift,  
                                                               int term_shift_num, 
                                                               int threshold, int top_num) {

  double min_mass = prec_minus_water_mass - max_shift - prec_error_tole;
  double max_mass = prec_minus_water_mass - min_shift + prec_error_tole;
  //LOG_DEBUG("min mass " << min_mass << " max mass " << max_mass);
  std::vector<int> row_begins = index_ptr->getProteoRowBegins();
  std::vector<int> row_ends = index_ptr->getProteoRowEnds();
  std::vector<double> proteo_minus_water_masses = index_ptr->getProteoMinusWaterMasses();
  std::vector<double> trunc_shifts = index_ptr->getTruncShifts();
  std::vector<int> rev_row_begins = rev_index_ptr->getProteoRowBegins();
  std::vector<int> rev_row_ends = rev_index_ptr->getProteoRowEnds();
  std::vector<double> rev_trunc_shifts = rev_index_ptr->getTruncShifts();

  ProtScorePtrVec proteo_scores;
  int bgn, end, rev_bgn, rev_end;
  int k_start, k_end;
  int best_score;
  double best_n_term_shift, best_c_term_shift;
  for (size_t i = 0; i < row_begins.size(); i++) {
    double proteo_mass = proteo_minus_water_masses[i];
    bgn = row_begins[i];
    end = row_ends[i];
    rev_bgn = rev_row_begins[i];
    rev_end = rev_row_ends[i];

    k_start = rev_end; 
    k_end = rev_end;
    best_score = 0;
    best_n_term_shift = 0;
    best_c_term_shift = 0;
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
      ProtScorePtr proteo_score = std::make_shared<ProtScore>(i, best_score, 
                                                              best_n_term_shift, 
                                                              best_c_term_shift);
      proteo_scores.push_back(proteo_score);
    }
  }

  ProtCandidatePtrVec proteins = ProtCandidate::geneResults(proteo_scores, threshold, top_num);
  // add trunc shifts
  addTruncShiftsWithRestriction(proteins, index_ptr, rev_index_ptr, 
                                proteo_minus_water_masses, 
                                scores, rev_scores, 
                                min_mass, max_mass, 
                                term_shift_num);
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
                                            int threshold, int top_num) {
  // if allowed mass shifts is small, use restricted search
  if ((max_shift - min_shift) <= 2000) {
    return findOneShiftTopProteinsWithRestrict(scores, rev_scores,
                                               index_ptr,rev_index_ptr,
                                               prec_minus_water_mass, prec_error_tole,
                                               min_shift, max_shift, 
                                               term_shift_num, threshold, top_num); 
  }
  else {
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
      cand_prots = findOneShiftTopProteinsWithRestrict(scores, rev_scores,
                                                       index_ptr,rev_index_ptr,
                                                       prec_minus_water_mass, prec_error_tole,
                                                       adjust_min_shift, adjust_max_shift, 
                                                       term_shift_num, threshold, top_num); 
    }
    ProtCandidatePtrVec unrestrict_cand_prots 
      =  findOneShiftTopProteinsWithoutRestrict(scores, rev_scores,
                                                index_ptr, rev_index_ptr,
                                                term_shift_num, threshold, top_num); 
    // concatenate
    cand_prots.insert(cand_prots.end(), unrestrict_cand_prots.begin(), unrestrict_cand_prots.end());
    return cand_prots;
  }
}

ProtCandidatePtrVec findVarPtmTopProteins(std::vector<short> &scores,
                                          std::vector<short> &rev_scores,
                                          MassMatchPtr index_ptr,
                                          MassMatchPtr rev_index_ptr,
                                          double prec_minus_water_mass,
                                          double prec_error_tole, 
                                          std::vector<double> &ptm_shifts,
                                          int threshold, int top_num) {

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
        double mass = proteo_minus_water_masses[i] + n_trunc_shift + c_trunc_shift; 
        if (mass < 0) {
          continue;
        }
        double score = scores[prec_match_rows[j]] + rev_scores[rev_prec_match_rows[k]];
        bool mass_match = false;
        for (size_t s = 0; s < ptm_shifts.size(); s++) {
          double shift_mass = mass + ptm_shifts[s];
          //LOG_DEBUG("shift mass " << shift_mass << " min mass " << min_mass << " max mass " << max_mass);
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
      ProtScorePtr proteo_score = std::make_shared<ProtScore>(i, best_score, 
                                                              best_n_trunc_shift, 
                                                              best_c_trunc_shift);
      proteo_scores.push_back(proteo_score);
    }
  }


  ProtCandidatePtrVec proteins = ProtCandidate::geneResults(proteo_scores, threshold, top_num);
  return proteins;
}

// Find top proteins for diagonal filtering
ProtCandidatePtrVec findDiagTopProteins(std::vector<short> &scores,
                                        std::vector<int> &proteo_row_begins,
                                        std::vector<int> &proteo_row_ends,
                                        int threshold, int top_num) {
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
  return ProtCandidate::geneResults(diag_scores, threshold, top_num);
}


} // namespace mass_match_util

} // namespace toppic 
