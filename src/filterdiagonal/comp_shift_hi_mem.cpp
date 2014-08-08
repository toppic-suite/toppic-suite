#include <iostream>

#include "base/logger.hpp"
#include "filterdiagonal/comp_shift_hi_mem.hpp"

namespace prot {

CompShiftHiMem::CompShiftHiMem(const ProteoformPtrVec &proteo_ptrs,
                               PtmFastFilterMngPtr mng_ptr){
  scale_ = mng_ptr->ptm_fast_filter_scale_;
  LOG_DEBUG("Scale: " + scale_);
  LOG_DEBUG("Proteoform point number: " + proteo_ptrs.size());

  col_num_ = mng_ptr->max_proteoform_mass * scale_;
  initProteoformBeginEnds(proteo_ptrs);
  initIndexes(proteo_ptrs);

  LOG_DEBUG("column number: " << col_num_);
  LOG_DEBUG("row number: " << row_num_);
  LOG_DEBUG("indexes size: "<< sizeof(col_indexes_)/sizeof(col_indexes_[0]));
}

CompShiftHiMem::~CompShiftHiMem(){
  delete proteo_row_begins_;
  delete row_proteo_ids_;
  delete col_index_begins_;
  delete col_index_ends_;
  delete col_indexes_;
}

void CompShiftHiMem::initProteoformBeginEnds(const ProteoformPtrVec &proteo_ptrs){
  proteo_row_begins_ = new int[proteo_ptrs.size()];
  int* proteo_row_ends;
  proteo_row_ends = new int[proteo_ptrs.size()];
  int  pnt = 0;
  for(size_t i=0; i< proteo_ptrs.size(); i++){
    proteo_row_begins_[i] = pnt;
    proteo_row_ends[i] = pnt + proteo_ptrs[i]->getResSeqPtr()->getLen() - 1;
    pnt += proteo_ptrs[i]->getResSeqPtr()->getLen();
  }
  row_num_ = pnt;
  row_proteo_ids_ = new int[row_num_];
  for(size_t i =0; i<proteo_ptrs.size(); i++){
    for(int j= proteo_row_begins_[i]; j<= proteo_row_ends[i];j++){
      row_proteo_ids_[j] = i;
    }
  }
  delete proteo_row_ends;
}

inline void CompShiftHiMem::updateColumnMatchNums(ProteoformPtr proteo_ptr, 
                                                  int* col_match_nums) {
  std::vector<int> masses = proteo_ptr->getBpSpecPtr()->getScaledPrmMasses(scale_);
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

void CompShiftHiMem::initIndexes(const ProteoformPtrVec &proteo_ptrs){
  int* col_match_nums = new int[col_num_];
  int* col_index_pnts = new int[col_num_];
  col_index_begins_ = new int[col_num_];
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

  col_indexes_ = new int[pnt];

  for(size_t i=0; i<proteo_ptrs.size(); i++){
    if (i/1000*1000 == i) {
      LOG_DEBUG("preprocessing proteoform " << i);
    }
    std::vector<int> masses  
        = proteo_ptrs[i]->getBpSpecPtr()->getScaledPrmMasses(scale_);
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
  delete col_match_nums;
  delete col_index_pnts;
}


std::vector<std::pair<int,int>> CompShiftHiMem::compConvolution(
    const std::vector<int> &masses,int bgn_pos,int num){

  short* scores = new short[row_num_];

  int begin_index;
  int end_index;
  int m;

  for (size_t i = bgn_pos+1; i<masses.size(); i++){
    m = masses[i]-masses[bgn_pos];
    if(m>=col_num_){
      break;
    }
    if(m>0){
      begin_index = col_index_begins_[m-1];
    }
    else{
      begin_index = col_index_begins_[m];
    }
    end_index = col_index_ends_[m+1];
    for(int j = begin_index;j<end_index;j++){
      scores[col_indexes_[j]]++;
    }
  }
  std::vector<std::pair<int,int>> results = getShiftScores(scores, num);
  delete scores;
  return results;
}

std::vector<std::pair<int,int>> CompShiftHiMem::compConvolution(
    const std::vector<int> &masses, const std::vector<int> &errors,int bgn_pos,int num){

  short* scores = new short[row_num_];

  int begin_index;
  int end_index;
  int m;

  for(size_t i =bgn_pos+1; i<masses.size(); i++){

    m = masses[i]-masses[bgn_pos];
    // m - errors[i] performs better than m - errors[i] -  errors[bgn_pos]
    int left = m-errors[i];
    if(left < 0){
      left=0;
    }
    int right = m+errors[i];
    if(right >= col_num_){
      break;
    }
    begin_index = col_index_begins_[left];
    end_index= col_index_ends_[right];
    for(int j=begin_index;j<=end_index;j++){
      scores[col_indexes_[j]]++;
    }
  }
  std::vector<std::pair<int,int>> results = getShiftScores(scores, num);
  delete scores;
  return results;
}

inline std::vector<std::pair<int,int>> CompShiftHiMem::getShiftScores(short* scores,int num){
  short* top_scores = new short[num];
  int* top_rows = new int[num];
  memset(top_rows, -1, num * sizeof(int));

  int last_scr = top_scores[num-1];
  for(int i=0; i<row_num_; i++){
    if(scores[i] <= last_scr){
      continue;
    }
    for(int j=num-2;j>=0;j--){
      if(scores[i]>top_scores[j]){
        top_scores[j+1] = top_scores[j];
        top_rows[j+1] = top_rows[j];
      }
      else{
        top_scores[j+1] = scores[i];
        top_rows[j+1] = i;
        break;
      }
    }
    if(scores[i]>top_scores[0]){
      top_scores[0] = scores[i];
      top_rows[0]=i;
    }
    last_scr = top_scores[num-1];
  }
  int output_num =0;
  for(int i=0;i<num;i++){
    if(top_rows[i]>=0){
      output_num++;
    }
  }
  std::vector<std::pair<int,int>> results;
  for(int i=0;i<output_num;i++){
    std::pair<int,int> proteo_score(row_proteo_ids_[top_rows[i]], top_scores[i]);
    results.push_back(proteo_score);
  }
  delete top_scores;
  delete top_rows;
  return results;
}


} /* namespace prot */
