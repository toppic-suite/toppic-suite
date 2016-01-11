#include <iostream>
#include <algorithm>

#include "base/logger.hpp"
#include "base/prot_mod.hpp"
#include "base/prot_mod_util.hpp"
#include "base/mass_constant.hpp"
#include "zeroptmfilter/comp_shift.hpp"

namespace prot {

CompShift::CompShift(const ProteoformPtrVec &proteo_ptrs, int scale,
                     double max_proteoform_mass, ProtModPtrVec prot_mod_ptr_vec,
                     bool use_rev) {
  scale_ = scale;
  LOG_DEBUG("Scale: " << scale_);
  LOG_DEBUG("Proteoform number: " << proteo_ptrs.size());

  col_num_ = max_proteoform_mass * scale_;
  proteo_num_ = proteo_ptrs.size();
  prot_mod_ptr_vec_ = prot_mod_ptr_vec;

  LOG_DEBUG("column number: " << col_num_);
  LOG_DEBUG("row number: " << row_num_);

  LOG_DEBUG("start init");
  initProteoformBeginEnds(proteo_ptrs);
  LOG_DEBUG("init indexes");
  initIndexes(proteo_ptrs);
  if (use_rev) {
    LOG_DEBUG("init rev indexes");
    initRevIndexes(proteo_ptrs);
  }

}

CompShift::~CompShift(){
}

inline void CompShift::initProteoformBeginEnds(const ProteoformPtrVec &proteo_ptrs){
  //no need to init
  LOG_DEBUG("proteome size " << proteo_ptrs.size());
  proteo_row_begins_.resize(proteo_ptrs.size());
  proteo_row_ends_.resize(proteo_ptrs.size());
  int  pnt = 0;
  LOG_DEBUG("start iteration");
  for(size_t i=0; i< proteo_ptrs.size(); i++){
    proteo_row_begins_[i] = pnt;
    int len = proteo_ptrs[i]->getResSeqPtr()->getLen() ;
    ResiduePtrVec res_ptr_vec = proteo_ptrs[i]->getResSeqPtr()->getResidues();
    ProtModPtr prot_mod_ptr = ProtModUtil::findNME_Acetylation(prot_mod_ptr_vec_, res_ptr_vec);
    acet_mods_.push_back(prot_mod_ptr);
    if (acet_mods_[i] != nullptr) {
      len++;
      //LOG_DEBUG("index " << i << " not nullptr");
    }
    proteo_row_ends_[i] = pnt + len - 1;
    pnt += len;
  }
  //LOG_DEBUG("end iteration");
  row_num_ = pnt;
  row_proteo_ids_.resize(row_num_);
  n_trunc_shifts_.resize(row_num_);
  c_trunc_shifts_.resize(row_num_);
  for(size_t i =0; i<proteo_ptrs.size(); i++){
    for(int j= proteo_row_begins_[i]; j<= proteo_row_ends_[i];j++){
      row_proteo_ids_[j] = i;
    }
    std::vector<int> masses = proteo_ptrs[i]->getBpSpecPtr()->getScaledPrmMasses(scale_);
    std::vector<double> double_masses = proteo_ptrs[i]->getBpSpecPtr()->getPrmMasses();
    if (acet_mods_[i] != nullptr) {
      int ace_mass = (int)std::round(- acet_mods_[i]->getProtShift() * scale_);
      masses.push_back(ace_mass);
      std::sort(masses.begin(), masses.end(),std::less<int>()); 
      double double_ace_mass = - acet_mods_[i]->getProtShift();
      double_masses.push_back(double_ace_mass);
      std::sort(double_masses.begin(), double_masses.end(),std::less<double>()); 
    }
    for(int j= proteo_row_begins_[i]; j<= proteo_row_ends_[i];j++){
      int pos = j - proteo_row_begins_[i];
      n_trunc_shifts_[j] = - double_masses[pos];
    }
    std::sort(double_masses.begin(), double_masses.end(),std::greater<double>()); 
    for(int j= proteo_row_begins_[i]; j<= proteo_row_ends_[i];j++){
      int pos = j - proteo_row_begins_[i];
      c_trunc_shifts_[j] = double_masses[pos] - double_masses[0];
      //LOG_DEBUG("c_trum shift " << c_trunc_shifts_[j]);
    }
  }
}

inline void CompShift::updateColumnMatchNums(ProteoformPtr proteo_ptr, ProtModPtr acet_mod, 
                                                    std::vector<int> &col_match_nums) {
  std::vector<int> masses = proteo_ptr->getBpSpecPtr()->getScaledPrmMasses(scale_);
  if (acet_mod != nullptr) {
    int ace_mass = (int)std::round(- acet_mod->getProtShift() * scale_);
    masses.push_back(ace_mass);
    std::sort(masses.begin(), masses.end(),std::less<int>()); 
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

inline void CompShift::initIndexes(const ProteoformPtrVec &proteo_ptrs){
  LOG_DEBUG("column num " << col_num_);
  std::vector<int> col_match_nums (col_num_, 0); 
  // no need to initalize 
  std::vector<int> col_index_pnts (col_num_);
  col_index_begins_.resize(col_num_);
  col_index_ends_.resize(col_num_);

  for(size_t i =0; i<proteo_ptrs.size(); i++){
    updateColumnMatchNums(proteo_ptrs[i], acet_mods_[i], col_match_nums);
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

  for(size_t i=0; i<proteo_ptrs.size(); i++){
    if (i/1000*1000 == i) {
      LOG_DEBUG("preprocessing proteoform " << i);
    }
    std::vector<int> masses  
        = proteo_ptrs[i]->getBpSpecPtr()->getScaledPrmMasses(scale_);

    if (acet_mods_[i] != nullptr) {
      int ace_mass = (int)std::round(- acet_mods_[i]->getProtShift() * scale_);
      masses.push_back(ace_mass);
      std::sort(masses.begin(), masses.end(),std::less<int>()); 
    }
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

inline void CompShift::updateRevColumnMatchNums(ProteoformPtr proteo_ptr, ProtModPtr acet_mod, 
                                                std::vector<int> &col_match_nums) {
  std::vector<int> masses = proteo_ptr->getBpSpecPtr()->getScaledPrmMasses(scale_);
  //LOG_DEBUG("mass lenth " << masses.size());
  if (acet_mod != nullptr) {
    int ace_mass = (int)std::round(- acet_mod->getProtShift() * scale_);
    masses.push_back(ace_mass);
  }
  std::sort(masses.begin(), masses.end(),std::greater<int>() ); 
  for (size_t bgn = 0; bgn < masses.size() - 1; bgn++) {
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


inline void CompShift::initRevIndexes(const ProteoformPtrVec &proteo_ptrs){
  LOG_DEBUG("start init rev col num " << col_num_);
  std::vector<int> rev_col_match_nums(col_num_, 0);
  // no need to initalize
  std::vector<int> rev_col_index_pnts(col_num_);
  rev_col_index_begins_.resize(col_num_);
  rev_col_index_ends_.resize(col_num_);

  for(size_t i =0; i<proteo_ptrs.size(); i++){
    updateRevColumnMatchNums(proteo_ptrs[i], acet_mods_[i], rev_col_match_nums);
  }
  LOG_DEBUG("update ended");

  int pnt = 0;
  for(int i=0; i< col_num_; i++){
    rev_col_index_begins_[i] = pnt;
    rev_col_index_pnts[i] = pnt;
    rev_col_index_ends_[i] = pnt + rev_col_match_nums[i]-1;
    pnt += rev_col_match_nums[i];
  }
  // no need to init
  rev_col_indexes_.resize(pnt);
  LOG_DEBUG("indexes size: "<< pnt);

  for(size_t i=0; i<proteo_ptrs.size(); i++){
    if (i/1000*1000 == i) {
      LOG_DEBUG("preprocessing proteoform " << i);
    }
    std::vector<int> masses  
        = proteo_ptrs[i]->getBpSpecPtr()->getScaledPrmMasses(scale_);

    if (acet_mods_[i] != nullptr) {
      int ace_mass = (int)std::round(- acet_mods_[i]->getProtShift() * scale_);
      masses.push_back(ace_mass);
    }
    std::sort(masses.begin(), masses.end(), std::greater<int>()); 
    for (size_t bgn=0; bgn < masses.size()-1; bgn++) {
      for (size_t cur = bgn+1; cur < masses.size(); cur++) {
        int diff = -(masses[cur] - masses[bgn]);
        
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
}

void CompShift::compScores(const std::vector<std::pair<int,int>> &pref_mass_errors,
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

void CompShift::compScores(const std::vector<std::pair<int,int>> &pref_mass_errors,
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

void CompShift::compRevScores(const std::vector<std::pair<int,int>> &suff_mass_errors,
                              std::vector<short> &rev_scores) {
  int begin_index;
  int end_index;
  int m;
  for(size_t i = 0; i < suff_mass_errors.size(); i++){
    m = suff_mass_errors[i].first;
    //LOG_DEBUG("REV_SP MASS " << m);
    int left = m-suff_mass_errors[i].second;
    //LOG_DEBUG("LEFT " << left);
    if(left < 0){
      left=0;
    }
    int right = m + suff_mass_errors[i].second;
    //LOG_DEBUG("RIGHT " << right);
    if (right < 0 || right >= col_num_) {
      continue;
    }
    begin_index = rev_col_index_begins_[left];
    end_index= rev_col_index_ends_[right];
    //LOG_DEBUG("begin index " << begin_index << " end index " << end_index);
    for(int j=begin_index;j<=end_index;j++){
      rev_scores[rev_col_indexes_[j]]++;
      //LOG_DEBUG("REV ROW INDEX " << rev_col_indexes_[j] << " rev score " << rev_scores[rev_col_indexes_[j]]);
    }
  }
}


/*
void CompShift::compRevScores(const std::vector<std::pair<int,int>> &mass_errors,
                              std::vector<short> &rev_scores) {
  int begin_index;
  int end_index;
  int m;
  for(size_t i = 0; i < mass_errors.size(); i++){
    m = mass_errors[i].first - MassConstant::getWaterMass() * scale_;
    //LOG_DEBUG("REV_SP MASS " << m);
    int left = m-mass_errors[i].second;
    //LOG_DEBUG("LEFT " << left);
    if(left < 0){
      left=0;
    }
    int right = m + mass_errors[i].second;
    //LOG_DEBUG("RIGHT " << right);
    if (right < 0 || right >= col_num_) {
      continue;
    }
    begin_index = rev_col_index_begins_[left];
    end_index= rev_col_index_ends_[right];
    //LOG_DEBUG("begin index " << begin_index << " end index " << end_index);
    for(int j=begin_index;j<=end_index;j++){
      rev_scores[rev_col_indexes_[j]]++;
      //LOG_DEBUG("REV ROW INDEX " << rev_col_indexes_[j] << " rev score " << rev_scores[rev_col_indexes_[j]]);
    }
  }
}
*/

void CompShift::compZeroPtmConvolution(const std::vector<std::pair<int,int>> &pref_mass_errors, 
                                       const std::vector<std::pair<int,int>> &suff_mass_errors,
                                       std::pair<int,int> &prec_minus_water_mass_error, 
                                       int comp_num, int pref_suff_num, int inte_num) {
  std::vector<short> scores(row_num_, 0);
  compScores(pref_mass_errors, scores);
  std::vector<short> rev_scores(row_num_, 0);
  compRevScores(suff_mass_errors, rev_scores);
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
  if(right >= 0 || right < col_num_){
    // update scores
    begin_index = col_index_begins_[left];
    end_index= col_index_ends_[right];
    for(int j=begin_index;j<=end_index;j++){
      scores[col_indexes_[j]] += PRECURSOR_MATCH_SCORE;
    }
    // update rev scores
    begin_index = rev_col_index_begins_[left];
    end_index= rev_col_index_ends_[right];
    //LOG_DEBUG("rev begin index " << begin_index << " rev end index " << end_index);
    for(int j=begin_index;j<=end_index;j++){
      rev_scores[rev_col_indexes_[j]] += PRECURSOR_MATCH_SCORE;
      //LOG_DEBUG("rev row index " << rev_col_indexes_[j] << " rev score " << rev_scores[rev_col_indexes_[j]]);
    }
  }
  double threshold = PRECURSOR_MATCH_SCORE * 2.0 + 4.0;
  findTopScores(scores, rev_scores, threshold, comp_num, pref_suff_num, inte_num, false, 0);
}

void CompShift::compOnePtmConvolution(const std::vector<std::pair<int,int>> &pref_mass_errors, 
                                      const std::vector<std::pair<int,int>> &suff_mass_errors,
                                      int comp_num, int pref_suff_num, int inte_num, int shift_num) {
  /*
  LOG_DEBUG("mass number " << mass_errors.size());
  for (size_t i = 0; i < mass_errors.size(); i++) {
    LOG_DEBUG(i << " mass " << mass_errors[i].first << " error " << mass_errors[i].second);
  }
  */
  std::vector<short> scores(row_num_, 0);
  compScores(pref_mass_errors, scores);
  std::vector<short> rev_scores(row_num_, 0);
  compRevScores(suff_mass_errors, rev_scores);
  double threshold = 4.0; 
  findTopScores(scores, rev_scores, threshold, comp_num, pref_suff_num, inte_num, true, shift_num);
}

void CompShift::compDiagConvolution(const std::vector<std::pair<int,int>> &mass_errors, 
                                    int start, int top_num) {
  std::vector<short> scores(row_num_, 0);
  compScores(mass_errors, start, scores);
  findTopDiagScores(scores, top_num);
}

inline bool cmpScore(const std::pair<int, int> &a, const std::pair<int, int> &b) {
    return a.second > b.second;
}


inline void addResults(FilterProteinPtrVec &results, 
                       std::vector<std::pair<int,int>> &single_type_results, 
                       double threshold, int single_type_num) {

  results.clear();
  std::sort(single_type_results.begin(), single_type_results.end(), cmpScore);
  int output_num =0;
  for(int i=0;i< single_type_num;i++){
    if (i >= (int)single_type_results.size()) {
      break;
    }
    //LOG_DEBUG("rank " << i << " score " << single_type_results[i].second);
    if(single_type_results[i].second >= threshold){
      output_num++;
    }
    else {
      break;
    }
  }
  for(int i=0;i<output_num;i++){
    int prot_id = single_type_results[i].first;
    int score = single_type_results[i].second;
    FilterProteinPtr prot_ptr = FilterProteinPtr(new FilterProtein(prot_id, score));
    results.push_back(prot_ptr);
  }
}

void CompShift::addNTruncShifts(FilterProteinPtrVec &prot_ptrs, 
                                std::vector<short> &scores, int shift_num) {
  for (size_t i = 0; i < prot_ptrs.size(); i++) {
    FilterProteinPtr prot_ptr = prot_ptrs[i];
    int prot_id = prot_ptr->getProteinId();
    int bgn = proteo_row_begins_[prot_id];
    int end = proteo_row_ends_[prot_id];
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
      shifts.push_back(n_trunc_shifts_[pos_scores[j].first]);
    }
    prot_ptr->setNTermShifts(shifts);
  }
}

void CompShift::addCTruncShifts(FilterProteinPtrVec &prot_ptrs, 
                                std::vector<short> &scores, int shift_num) {
  for (size_t i = 0; i < prot_ptrs.size(); i++) {
    FilterProteinPtr prot_ptr = prot_ptrs[i];
    int prot_id = prot_ptr->getProteinId();
    int bgn = proteo_row_begins_[prot_id];
    int end = proteo_row_ends_[prot_id];
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
      shifts.push_back(c_trunc_shifts_[pos_scores[j].first]);
    }
    prot_ptr->setCTermShifts(shifts);
  }
}

inline void CompShift::findTopScores(std::vector<short> &scores, std::vector<short> &rev_scores, 
                                     double threshold, int comp_num, int pref_suff_num, int inte_num,
                                     bool add_shifts, int shift_num) {
  std::vector<std::pair<int,int>> comp_proteo_scores;
  std::vector<std::pair<int,int>> pref_proteo_scores;
  std::vector<std::pair<int,int>> suff_proteo_scores;
  std::vector<std::pair<int,int>> internal_proteo_scores;
  for (int i = 0; i < proteo_num_; i++) {
    int bgn = proteo_row_begins_[i];
    int end = proteo_row_ends_[i];
    int pref = bgn + 1;
    if (acet_mods_[i]) {
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
    //LOG_DEBUG("begin " << bgn << " end " << end << " rev 0 " << rev_scores[0]);
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
      if (j <= suff && rev_scores[j] > best_suff_score) {
        best_suff_score = rev_scores[j];
      }
    }
    //LOG_DEBUG("best pref score "  << best_pref_score << " best suffix score " << best_suff_score);
    std::pair<int,int> comp_proteo_score(i, best_pref_score + best_suff_score);
    comp_proteo_scores.push_back(comp_proteo_score);
    std::pair<int,int> pref_proteo_score(i, best_pref_score + best_rev_score);
    pref_proteo_scores.push_back(pref_proteo_score);
    std::pair<int,int> suff_proteo_score(i, best_score + best_suff_score);
    suff_proteo_scores.push_back(suff_proteo_score);
    std::pair<int,int> internal_proteo_score(i, best_score + best_rev_score);
    internal_proteo_scores.push_back(internal_proteo_score);
  }
  //LOG_DEBUG("num " << num << " Single type num " << single_type_num);
  addResults(top_comp_prots_, comp_proteo_scores, threshold, comp_num);
  addResults(top_pref_prots_, pref_proteo_scores, threshold, pref_suff_num);
  addResults(top_suff_prots_, suff_proteo_scores, threshold, pref_suff_num);
  addResults(top_internal_prots_, internal_proteo_scores, threshold, inte_num);
  if (add_shifts) {
    addNTruncShifts(top_suff_prots_, scores, shift_num);
    addNTruncShifts(top_internal_prots_, scores, shift_num);
    addCTruncShifts(top_pref_prots_, rev_scores, shift_num);
    addCTruncShifts(top_internal_prots_, rev_scores, shift_num);
  }
}

inline void CompShift::findTopDiagScores(std::vector<short> &scores, int num) {
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
  addResults(top_diag_prots_, diag_scores, threshold, num);
}

} /* namespace prot */
