#include <algorithm>

#include "base/logger.hpp"
#include "zeroptmfilter/mass_match_util.hpp"

namespace prot {

FilterProteinPtrVec MassMatchUtil::findTopProteins(std::vector<short> &scores, 
                                                   std::vector<int> &proteo_row_begins,
                                                   std::vector<int> &proteo_row_ends,
                                                   int threshold, int num) {
  std::vector<std::pair<int,int>> diag_scores;
  for (size_t i = 0; i < proteo_row_begins.size(); i++) {
    int bgn = proteo_row_begins[i];
    int end = proteo_row_ends[i];
    int best_score = 0;
    //LOG_DEBUG("begin " << bgn << " end " << end << " rev 0 " << rev_scores[0]);
    for (int j = bgn; j <= end; j++) {
      if (scores[j] > best_score) {
        best_score = scores[j];
      }
    }
    std::pair<int,int> diag_score(i, best_score);
    diag_scores.push_back(diag_score);
  }
  //LOG_DEBUG("num " << num << " Single type num " << single_type_num);
  return FilterProtein::geneResults(diag_scores, threshold, num);
  //LOG_DEBUG("top diag size " << top_diag_prots_.size());
}

inline bool cmpScore(const std::pair<int, int> &a, const std::pair<int, int> &b) {
    return a.second > b.second;
}

void addTruncShifts(FilterProteinPtrVec &prot_ptrs, MassMatchPtr index_ptr, 
                    std::vector<short> &scores, int shift_num, bool n_term) {
  std::vector<int> row_begins = index_ptr->getProteoRowBegins(); 
  std::vector<int> row_ends = index_ptr->getProteoRowEnds();
  std::vector<double> trunc_shifts = index_ptr->getTruncShifts();
  for (size_t i = 0; i < prot_ptrs.size(); i++) {
    FilterProteinPtr prot_ptr = prot_ptrs[i];
    int prot_id = prot_ptr->getProteinId();
    int bgn = row_begins[prot_id];
    int end = row_ends[prot_id];
    std::vector<std::pair<int, int>> pos_scores;
    pos_scores.resize(end-bgn + 1);
    for (int j = bgn; j <= end; j++) {
      std::pair<int, int> cur_pos_score (j, scores[j]);
      pos_scores[j-bgn] = cur_pos_score;
    }
    std::sort(pos_scores.begin(), pos_scores.end(), cmpScore);
    std::vector<double> shifts;
    for (size_t j = 0; j < pos_scores.size(); j++) {
      if ((int)j >= shift_num) {
        break;
      }
      shifts.push_back(trunc_shifts[pos_scores[j].first]);
    }
    if (n_term) {
      prot_ptr->setNTermShifts(shifts);
    }
    else {
      prot_ptr->setCTermShifts(shifts);
    }
  }
}

FilterProteinPtrVec MassMatchUtil::findTopProteins(std::vector<short> &scores, 
                                                   std::vector<short> &rev_scores, 
                                                   MassMatchPtr index_ptr,
                                                   MassMatchPtr rev_index_ptr,
                                                   double threshold, int num,
                                                   bool add_shifts, int shift_num) {
  std::vector<int> row_begins = index_ptr->getProteoRowBegins(); 
  std::vector<int> row_ends = index_ptr->getProteoRowEnds();
  std::vector<int> rev_row_begins = rev_index_ptr->getProteoRowBegins(); 
  std::vector<int> rev_row_ends = rev_index_ptr->getProteoRowEnds();
  std::vector<std::pair<int,int>> proteo_scores;
  for (size_t i = 0; i < row_begins.size(); i++) {
    int bgn = row_begins[i];
    int end = row_ends[i];
    int best_score = 0;
    for (int j = bgn; j <= end; j++) {
      //LOG_DEBUG("begin " << bgn << " end " << end << " score " << scores[j]);
      if (scores[j] > best_score) {
        best_score = scores[j];
      }
    }
   
    int rev_bgn = rev_row_begins[i];
    int rev_end = rev_row_ends[i];
    int rev_best_score = 0;
    for (int j = rev_bgn; j <= rev_end; j++) {
      //LOG_DEBUG("begin " << bgn << " end " << end << " rev_score " << rev_scores[j]);
      if (rev_scores[j] > rev_best_score) {
        rev_best_score = rev_scores[j];
      }
    }
    //LOG_DEBUG("best score " << best_score << " rev best score " << rev_best_score);

    std::pair<int,int> proteo_score(i, best_score + rev_best_score);
    proteo_scores.push_back(proteo_score);
  }

  FilterProteinPtrVec proteins = FilterProtein::geneResults(proteo_scores, threshold, num);
  if (add_shifts) {
    bool n_term = true;
    addTruncShifts(proteins, index_ptr, scores, shift_num, n_term);
    n_term = false;
    addTruncShifts(proteins, rev_index_ptr, rev_scores, shift_num, n_term);
  }
  return proteins;
}

} /* namespace prot */
