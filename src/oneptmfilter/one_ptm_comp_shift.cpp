#include <iostream>
#include <algorithm>

#include "base/logger.hpp"
#include "base/prot_mod.hpp"
#include "base/mass_constant.hpp"
#include "oneptmfilter/one_ptm_comp_shift.hpp"

namespace prot {

OnePtmCompShift::OnePtmCompShift(const ProteoformPtrVec &proteo_ptrs,
                                 OnePtmFilterMngPtr mng_ptr) {
  scale_ = mng_ptr->ptm_fast_filter_scale_;
  LOG_DEBUG("Scale: " << scale_);
  LOG_DEBUG("Proteoform number: " << proteo_ptrs.size());

  col_num_ = mng_ptr->max_proteoform_mass * scale_;
  proteo_num_ = proteo_ptrs.size();
  acetylation_ = containNME_ACETYLATION(mng_ptr->prsm_para_ptr_->getAllowProtModPtrVec());
  initProteoformBeginEnds(proteo_ptrs);
  initIndexes(proteo_ptrs);
  initRevIndexes(proteo_ptrs);

  LOG_DEBUG("column number: " << col_num_);
  LOG_DEBUG("row number: " << row_num_);
}

OnePtmCompShift::~OnePtmCompShift(){
  delete[] proteo_row_begins_;
  delete[] proteo_row_ends_;
  delete[] row_proteo_ids_;
  delete[] col_index_begins_;
  delete[] col_index_ends_;
  delete[] col_indexes_;
  delete[] rev_col_index_begins_;
  delete[] rev_col_index_ends_;
  delete[] rev_col_indexes_;
}

inline void OnePtmCompShift::initProteoformBeginEnds(const ProteoformPtrVec &proteo_ptrs){
  //no need to init
  proteo_row_begins_ = new int[proteo_ptrs.size()];
  //no need to init
  proteo_row_ends_ = new int[proteo_ptrs.size()];
  int  pnt = 0;
  for(size_t i=0; i< proteo_ptrs.size(); i++){
    proteo_row_begins_[i] = pnt;
    int len = proteo_ptrs[i]->getResSeqPtr()->getLen() ;
    if (acetylation_) {
      len++;
    }
    proteo_row_ends_[i] = pnt + len - 1;
    pnt += len;
  }
  row_num_ = pnt;
  //no need to init
  row_proteo_ids_ = new int[row_num_];
  for(size_t i =0; i<proteo_ptrs.size(); i++){
    for(int j= proteo_row_begins_[i]; j<= proteo_row_ends_[i];j++){
      row_proteo_ids_[j] = i;
    }
  }
}

inline void OnePtmCompShift::updateColumnMatchNums(ProteoformPtr proteo_ptr, 
                                                  int* col_match_nums) {
  std::vector<int> masses = proteo_ptr->getBpSpecPtr()->getScaledPrmMasses(scale_);
  if (acetylation_) {
    double ace_mass = - ProtModFactory::getProtModPtr_NME_ACETYLATION()->getProtShift();
    masses.push_back(ace_mass);
    std::sort(masses.begin(), masses.end(),std::less<double>()); 
  }
  for (size_t bgn = 0; bgn < masses.size(); bgn++) {
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

inline void OnePtmCompShift::initIndexes(const ProteoformPtrVec &proteo_ptrs){
  int* col_match_nums = new int[col_num_];
  // init col_match_nums
  memset(col_match_nums, 0, col_num_ * sizeof(int));
  // no need to init
  int* col_index_pnts = new int[col_num_];
  // no need to init
  col_index_begins_ = new int[col_num_];
  // no need to init
  col_index_ends_ = new int[col_num_];

  for(size_t i =0; i<proteo_ptrs.size(); i++){
    updateColumnMatchNums(proteo_ptrs[i], col_match_nums);
  }

  int pnt = 0;
  for(int i=0; i< col_num_; i++){
    col_index_begins_[i] = pnt;
    col_index_pnts[i] = pnt;
    col_index_ends_[i] = pnt + col_match_nums[i]-1;
    pnt += col_match_nums[i];
  }
  // no need to init
  col_indexes_ = new int[pnt];
  LOG_DEBUG("indexes size: "<< pnt);

  for(size_t i=0; i<proteo_ptrs.size(); i++){
    if (i/1000*1000 == i) {
      LOG_DEBUG("preprocessing proteoform " << i);
    }
    std::vector<int> masses  
        = proteo_ptrs[i]->getBpSpecPtr()->getScaledPrmMasses(scale_);
    if (acetylation_) {
      double ace_mass = - ProtModFactory::getProtModPtr_NME_ACETYLATION()->getProtShift();
      masses.push_back(ace_mass);
      std::sort(masses.begin(), masses.end(),std::less<double>()); 
    }
    for (size_t bgn=0; bgn < masses.size(); bgn++) {
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
  delete[] col_match_nums;
  delete[] col_index_pnts;
}

inline void OnePtmCompShift::updateRevColumnMatchNums(ProteoformPtr proteo_ptr, 
                                                      int* col_match_nums) {
  std::vector<int> masses = proteo_ptr->getBpSpecPtr()->getScaledPrmMasses(scale_);
  if (acetylation_) {
    double ace_mass = - ProtModFactory::getProtModPtr_NME_ACETYLATION()->getProtShift();
    masses.push_back(ace_mass);
  }
  std::sort(masses.begin(), masses.end(),std::greater<double>() ); 
  for (size_t bgn = 0; bgn < masses.size(); bgn++) {
    for (size_t cur = bgn + 1; cur < masses.size(); cur++) {
      int diff = -(masses[cur] - masses[bgn]);
      if (diff < col_num_) {
        col_match_nums[diff]++;
      }
      else {
        break;
      }
    }
  }
}


inline void OnePtmCompShift::initRevIndexes(const ProteoformPtrVec &proteo_ptrs){
  int* rev_col_match_nums = new int[col_num_];
  // init col_match_nums
  memset(rev_col_match_nums, 0, col_num_ * sizeof(int));
  // no need to init
  int* rev_col_index_pnts = new int[col_num_];
  // no need to init
  rev_col_index_begins_ = new int[col_num_];
  // no need to init
  rev_col_index_ends_ = new int[col_num_];

  for(size_t i =0; i<proteo_ptrs.size(); i++){
    updateRevColumnMatchNums(proteo_ptrs[i], rev_col_match_nums);
  }

  int pnt = 0;
  for(int i=0; i< col_num_; i++){
    rev_col_index_begins_[i] = pnt;
    rev_col_index_pnts[i] = pnt;
    rev_col_index_ends_[i] = pnt + rev_col_match_nums[i]-1;
    pnt += rev_col_match_nums[i];
  }
  // no need to init
  rev_col_indexes_ = new int[pnt];
  LOG_DEBUG("indexes size: "<< pnt);

  for(size_t i=0; i<proteo_ptrs.size(); i++){
    if (i/1000*1000 == i) {
      LOG_DEBUG("preprocessing proteoform " << i);
    }
    std::vector<int> masses  
        = proteo_ptrs[i]->getBpSpecPtr()->getScaledPrmMasses(scale_);
    if (acetylation_) {
      double ace_mass = - ProtModFactory::getProtModPtr_NME_ACETYLATION()->getProtShift();
      masses.push_back(ace_mass);
    }

    std::sort(masses.begin(), masses.end(), std::greater<double>()); 
    for (size_t bgn=0; bgn < masses.size(); bgn++) {
      for (size_t cur = bgn+1; cur < masses.size(); cur++) {
        int diff = -(masses[cur] - masses[bgn]);
        
        /*
        if (bgn == 0) {
          LOG_DEBUG("INT MASS " << diff);
        }
        */

        if (diff < col_num_) {
          rev_col_indexes_[rev_col_index_pnts[diff]] = proteo_row_begins_[i] + bgn;
          rev_col_index_pnts[diff]++;
        }
        else {
          break;
        }
      }
    }
  }
  delete[] rev_col_match_nums;
  delete[] rev_col_index_pnts;
}


std::vector<std::pair<int,int>> OnePtmCompShift::compConvolution(
    std::vector<int> &masses, std::vector<int> &errors, int num){

  short* scores = new short[row_num_];
  memset(scores, 0, row_num_ * sizeof(short));

  int begin_index;
  int end_index;
  int m;

  for(size_t i = 0; i<masses.size(); i++){

    m = masses[i];
    // m - errors[i] performs better than m - errors[i] -  errors[bgn_pos]
    int left = m-errors[i];
    if(left < 0){
      left=0;
    }
    int right = m+errors[i];
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

  short* rev_scores = new short[row_num_];
  memset(rev_scores, 0, row_num_ * sizeof(short));

  for(size_t i = 0; i<masses.size(); i++){

    m = masses[i] - MassConstant::getWaterMass() * scale_;
    //LOG_DEBUG("REV_SP MASS " << m);
    int left = m-errors[i];
    //LOG_DEBUG("LEFT " << left);
    if(left < 0){
      left=0;
    }
    int right = m + errors[i];
    //LOG_DEBUG("RIGHT " << right);
    if (right < 0 || right >= col_num_) {
      continue;
    }
    begin_index = rev_col_index_begins_[left];
    end_index= rev_col_index_ends_[right];
    for(int j=begin_index;j<=end_index;j++){
      rev_scores[rev_col_indexes_[j]]++;
      //LOG_DEBUG("REV ROW INDEX " << rev_col_indexes_[j] << " rev score " << rev_scores[rev_col_indexes_[j]]);
    }
  }

  std::vector<std::pair<int,int>> results = getShiftScores(scores, rev_scores, num);
  delete[] scores;
  delete[] rev_scores;
  return results;
}

inline bool scoreCompare(const std::pair<int, int> &a, const std::pair<int, int> &b) {
    return a.second > b.second;
}

inline void addResults(std::vector<std::pair<int,int>> &results, std::vector<std::pair<int,int>> &single_type_results, 
                       int single_type_num) {

  std::sort(single_type_results.begin(), single_type_results.end(), scoreCompare);
  int output_num =0;
  for(int i=0;i< single_type_num;i++){
    if (i >= (int)single_type_results.size()) {
      break;
    }
    //LOG_DEBUG("rank " << i << " score " << single_type_results[i].second);
    if(single_type_results[i].second > 4){
      output_num++;
    }
  }
  for(int i=0;i<output_num;i++){
    results.push_back(single_type_results[i]);
  }
}

inline std::vector<std::pair<int,int>> OnePtmCompShift::getShiftScores(short* scores, short* rev_scores, int num){
  std::vector<std::pair<int,int>> comp_proteo_scores;
  std::vector<std::pair<int,int>> pref_proteo_scores;
  std::vector<std::pair<int,int>> suff_proteo_scores;
  std::vector<std::pair<int,int>> internal_proteo_scores;
  for (int i = 0; i < proteo_num_; i++) {
    int bgn = proteo_row_begins_[i];
    int end = proteo_row_ends_[i];
    int pref = bgn + 1;
    if (acetylation_) {
      pref = bgn + 2;
    }
    if (pref > end) {
      pref = end;
    }
    int suff = bgn;
    int best_score = 0;
    int best_pref_score = 0;
    int best_rev_score = 0;
    int best_suff_score = 0;
    for (int j = bgn; j <= end; j++) {
      if (scores[j] > best_score) {
        best_score = scores[j];
      }
      if (j <= pref && scores[j] > best_pref_score) {
        best_pref_score = scores[j];
      }
      if (rev_scores[j] > best_rev_score) {
        best_rev_score = rev_scores[j];
      }
      if (j <= suff && scores[j] > best_suff_score) {
        best_suff_score = rev_scores[j];
      }
    }
    std::pair<int,int> comp_proteo_score(i, best_pref_score + best_suff_score);
    comp_proteo_scores.push_back(comp_proteo_score);
    std::pair<int,int> pref_proteo_score(i, best_pref_score + best_rev_score);
    pref_proteo_scores.push_back(pref_proteo_score);
    std::pair<int,int> suff_proteo_score(i, best_score + best_suff_score);
    suff_proteo_scores.push_back(suff_proteo_score);
    std::pair<int,int> internal_proteo_score(i, best_score + best_rev_score);
    internal_proteo_scores.push_back(suff_proteo_score);
  }
  int single_type_num = num / 4;
  //LOG_DEBUG("num " << num << " Single type num " << single_type_num);
  std::vector<std::pair<int,int>> results;
  addResults(results, comp_proteo_scores, single_type_num);
  addResults(results, pref_proteo_scores, single_type_num);
  addResults(results, suff_proteo_scores, single_type_num);
  addResults(results, internal_proteo_scores, single_type_num);

  return results;
}


} /* namespace prot */
