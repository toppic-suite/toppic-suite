#include <iostream>
#include <algorithm>

#include "base/logger.hpp"
#include "base/prot_mod.hpp"
#include "base/prot_mod_util.hpp"
#include "base/mass_constant.hpp"
#include "zeroptmfilter/filter_protein.hpp"
#include "zeroptmfilter/mass_match.hpp"

namespace prot {

MassMatch::MassMatch(const ProteoformPtrVec &proteo_ptrs, int scale, 
                     double max_proteoform_mass, bool rev) {
  scale_ = scale;
  LOG_DEBUG("Scale: " << scale_);
  LOG_DEBUG("Proteoform number: " << proteo_ptrs.size());

  col_num_ = max_proteoform_mass * scale_;
  proteo_num_ = proteo_ptrs.size();
  LOG_DEBUG("column number: " << col_num_);
  LOG_DEBUG("row number: " << row_num_);

  std::vector<int> index_nums;
  for (size_t i = 0; i < proteo_ptrs.size(); i++) {
    index_nums.push_back(proteo_ptrs[i]->getResSeqPtr()->getLen());
  }

  LOG_DEBUG("start init");
  initProteoformBeginEnds(proteo_ptrs, index_nums);
  LOG_DEBUG("init indexes");
  initIndexes(proteo_ptrs, rev);
}

MassMatch::MassMatch(const ProteoformPtrVec &proteo_ptrs, 
                     const std::vector<std::vector<int>> &shift_2d,
                     int scale, double max_proteoform_mass, bool rev) {
  scale_ = scale;
  LOG_DEBUG("Scale: " << scale_);
  LOG_DEBUG("Proteoform number: " << proteo_ptrs.size());

  col_num_ = max_proteoform_mass * scale_;
  proteo_num_ = proteo_ptrs.size();

  LOG_DEBUG("column number: " << col_num_);
  LOG_DEBUG("row number: " << row_num_);
  std::vector<int> index_nums;
  for (size_t i = 0; i < proteo_ptrs.size(); i++) {
    index_nums.push_back(shift_2d[i].size());
  }

  LOG_DEBUG("start init");
  initProteoformBeginEnds(proteo_ptrs, index_nums);
  LOG_DEBUG("init indexes");
  initIndexes(proteo_ptrs, shift_2d, rev);
}

MassMatch::~MassMatch(){
}

inline void MassMatch::initProteoformBeginEnds(const ProteoformPtrVec &proteo_ptrs,
                                               const std::vector<int> &index_nums){
  //no need to init
  LOG_DEBUG("proteome size " << proteo_ptrs.size());
  proteo_row_begins_.resize(proteo_ptrs.size());
  proteo_row_ends_.resize(proteo_ptrs.size());
  int  pnt = 0;
  LOG_DEBUG("start iteration");
  for(size_t i=0; i< proteo_ptrs.size(); i++){
    proteo_row_begins_[i] = pnt;
    int len = index_nums[i];
    proteo_row_ends_[i] = pnt + len - 1;
    pnt += len;
  }
  //LOG_DEBUG("end iteration");
  row_num_ = pnt;
  row_proteo_ids_.resize(row_num_);
  for(size_t i =0; i<proteo_ptrs.size(); i++){
    for(int j= proteo_row_begins_[i]; j<= proteo_row_ends_[i];j++){
      row_proteo_ids_[j] = i;
    }
  }
}

inline void MassMatch::updateColumnMatchNums(ProteoformPtr proteo_ptr,
                                             std::vector<int> &col_match_nums,
                                             bool rev) {
  std::vector<int> masses;
  if (!rev) {
    masses = proteo_ptr->getBpSpecPtr()->getScaledPrmMasses(scale_);
  }
  else {
    masses = proteo_ptr->getBpSpecPtr()->getScaledSrmMasses(scale_);
  }

  for (size_t bgn = 0; bgn < masses.size() - 1; bgn++) {
    for (size_t cur = bgn + 1; cur < masses.size(); cur++) {
      int diff = masses[cur] - masses[bgn];
      if (diff < col_num_) {
        col_match_nums[diff]++;
      }
      else {
        break;
      }
    }
  }
}

inline void MassMatch::initIndexes(const ProteoformPtrVec &proteo_ptrs, bool rev){
  LOG_DEBUG("column num " << col_num_);
  std::vector<int> col_match_nums (col_num_, 0); 
  // no need to initalize 
  std::vector<int> col_index_pnts (col_num_);
  col_index_begins_.resize(col_num_);
  col_index_ends_.resize(col_num_);

  for(size_t i =0; i<proteo_ptrs.size(); i++){
    updateColumnMatchNums(proteo_ptrs[i], col_match_nums, rev);
  }

  int pnt = 0;
  for(int i=0; i< col_num_; i++){
    col_index_begins_[i] = pnt;
    col_index_pnts[i] = pnt;
    col_index_ends_[i] = pnt + col_match_nums[i]-1;
    pnt += col_match_nums[i];
  }
  // no need to init
  col_indexes_.resize(pnt, 0);
  LOG_DEBUG("indexes size: "<< pnt);
  trunc_shifts_.resize(row_num_);

  for(size_t i=0; i<proteo_ptrs.size(); i++){
    if (i/1000*1000 == i) {
      LOG_DEBUG("preprocessing proteoform " << i);
    }
    std::vector<int> masses;  
    std::vector<double> double_masses;
    if (!rev) {
      masses = proteo_ptrs[i]->getBpSpecPtr()->getScaledPrmMasses(scale_);
      double_masses = proteo_ptrs[i]->getBpSpecPtr()->getPrmMasses();
    }
    else {
      masses = proteo_ptrs[i]->getBpSpecPtr()->getScaledSrmMasses(scale_);
      double_masses = proteo_ptrs[i]->getBpSpecPtr()->getSrmMasses();
    }

    // init trunc_shifts
    for(int j= proteo_row_begins_[i]; j<= proteo_row_ends_[i];j++){
      int pos = j - proteo_row_begins_[i];
      trunc_shifts_[j] = - double_masses[pos];
    }

    // init col indexes
    for (size_t bgn=0; bgn < masses.size() - 1; bgn++) {
      for (size_t cur = bgn+1; cur < masses.size(); cur++) {
        int diff = masses[cur] - masses[bgn];
        if (diff < col_num_) {
          col_indexes_[col_index_pnts[diff]] = proteo_row_begins_[i] + bgn;
          col_index_pnts[diff]++;
        }
        else {
          break;
        }
      }
    }
  }
}

inline void MassMatch::updateColumnMatchNums(ProteoformPtr proteo_ptr,
                                             const std::vector<int> &int_shifts,
                                             std::vector<int> &col_match_nums, 
                                             bool rev) {
  std::vector<int> masses;
  if (!rev) {
    masses = proteo_ptr->getBpSpecPtr()->getScaledPrmMasses(scale_);
  }
  else {
    masses = proteo_ptr->getBpSpecPtr()->getScaledSrmMasses(scale_);
  }
  for (size_t s = 0; s < int_shifts.size(); s++)  {
    for (size_t cur = 1; cur < masses.size(); cur++) {
      int shift_mass = masses[cur] + int_shifts[s];
      if (shift_mass > 0 && shift_mass < col_num_) {
        col_match_nums[shift_mass]++;
      }
      else {
        break;
      }
    }
  }
}

inline void MassMatch::initIndexes(const ProteoformPtrVec &proteo_ptrs, 
                                   const std::vector<std::vector<int>> &shift_2d,
                                   bool rev){
  LOG_DEBUG("column num " << col_num_);
  std::vector<int> col_match_nums (col_num_, 0); 
  // no need to initalize 
  std::vector<int> col_index_pnts (col_num_);
  col_index_begins_.resize(col_num_);
  col_index_ends_.resize(col_num_);

  for(size_t i =0; i<proteo_ptrs.size(); i++){
    updateColumnMatchNums(proteo_ptrs[i], shift_2d[i], col_match_nums, rev);
  }

  int pnt = 0;
  for(int i=0; i< col_num_; i++){
    col_index_begins_[i] = pnt;
    col_index_pnts[i] = pnt;
    col_index_ends_[i] = pnt + col_match_nums[i]-1;
    pnt += col_match_nums[i];
  }
  // no need to init
  col_indexes_.resize(pnt, 0);
  LOG_DEBUG("indexes size: "<< pnt);
  trunc_shifts_.resize(row_num_);

  for(size_t i=0; i<proteo_ptrs.size(); i++){
    if (i/1000*1000 == i) {
      LOG_DEBUG("preprocessing proteoform " << i);
    }
    std::vector<int> masses;
    if (!rev) {
      masses = proteo_ptrs[i]->getBpSpecPtr()->getScaledPrmMasses(scale_);
    }
    else {
      masses = proteo_ptrs[i]->getBpSpecPtr()->getScaledSrmMasses(scale_);
    }

    for(int j= proteo_row_begins_[i]; j<= proteo_row_ends_[i];j++){
      int pos = j - proteo_row_begins_[i];
      trunc_shifts_[j] = shift_2d[i][pos];
    }

    for (size_t s=0; s < shift_2d[i].size(); s++) {
      for (size_t cur = 1; cur < masses.size(); cur++) {
        int shift_mass = masses[cur] + shift_2d[i][s];
        if (shift_mass > 0 && shift_mass < col_num_) {
          col_indexes_[col_index_pnts[shift_mass]] = proteo_row_begins_[i] + s;
          col_index_pnts[shift_mass]++;
        }
        else {
          break;
        }
      }
    }
  }
}

void MassMatch::compScores(const std::vector<std::pair<int,int>> &pref_mass_errors,
                           std::vector<short> &scores) {
  int begin_index;
  int end_index;
  int m;
  for(size_t i = 0; i<pref_mass_errors.size(); i++){
    m = pref_mass_errors[i].first;
    // m - errors[i] performs better than m - errors[i] -  errors[bgn_pos]
    int left = m - pref_mass_errors[i].second;
    if(left < 0){
      left=0;
    }
    int right = m + pref_mass_errors[i].second;
    if(right < 0 || right >= col_num_){
      continue;
    }
    //LOG_DEBUG("SP MASS " << m);
    begin_index = col_index_begins_[left];
    end_index= col_index_ends_[right];
    for(int j=begin_index;j<=end_index;j++){
      scores[col_indexes_[j]]++;
      //LOG_DEBUG("ROW INDEX " << col_indexes_[j] << " score " << scores[col_indexes_[j]]);
    }
  }
}

void MassMatch::compScores(const std::vector<std::pair<int,int>> &pref_mass_errors,
                           int start, std::vector<short> &scores) {
  int begin_index;
  int end_index;
  int m;
  for(size_t i = start; i<pref_mass_errors.size(); i++){
    m = pref_mass_errors[i].first - pref_mass_errors[start].first;
    // m - errors[i] performs better than m - errors[i] -  errors[bgn_pos]
    int left = m - pref_mass_errors[i].second;
    if(left < 0){
      left=0;
    }
    int right = m + pref_mass_errors[i].second;
    if(right < 0 || right >= col_num_){
      continue;
    }
    //LOG_DEBUG("SP MASS " << m);
    begin_index = col_index_begins_[left];
    end_index= col_index_ends_[right];
    for(int j=begin_index;j<=end_index;j++){
      scores[col_indexes_[j]]++;
      //LOG_DEBUG("ROW INDEX " << col_indexes_[j] << " score " << scores[col_indexes_[j]]);
    }
  }
}

void MassMatch::compMatchScores(const std::vector<std::pair<int,int>> &pref_mass_errors,
                                const std::pair<int,int> &prec_minus_water_mass_error, 
                                std::vector<short> &scores) {
  compScores(pref_mass_errors, 0, scores);
  // precursor mass 
  int begin_index, end_index;
  int m = prec_minus_water_mass_error.first;
  // m - errors[i] performs better than m - errors[i] -  errors[bgn_pos]
  int left = m - prec_minus_water_mass_error.second;
  if(left < 0){
    left=0;
  }
  int right = m + prec_minus_water_mass_error.second;
  //LOG_DEBUG("prec left " << left << " pref right " << right);
  if(right >= 0 && right < col_num_){
    // update scores
    begin_index = col_index_begins_[left];
    end_index= col_index_ends_[right];
    for(int j=begin_index;j<=end_index;j++){
      scores[col_indexes_[j]] += getPrecursorMatchScore();
    }
  }
}

/*
FilterProteinPtrVec MassMatch::findTopProteins(std::vector<short> &scores, int num) {
  std::vector<std::pair<int,int>> diag_scores;
  for (int i = 0; i < proteo_num_; i++) {
    int bgn = proteo_row_begins_[i];
    int end = proteo_row_ends_[i];
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
  double threshold = 4.0;
  return FilterProtein::geneResults(diag_scores, threshold, num);
  //LOG_DEBUG("top diag size " << top_diag_prots_.size());
}
*/

} /* namespace prot */
