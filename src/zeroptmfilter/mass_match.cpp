#include <iostream>
#include <algorithm>

#include "base/logger.hpp"
#include "base/prot_mod.hpp"
#include "base/prot_mod_util.hpp"
#include "base/mass_constant.hpp"
#include "zeroptmfilter/filter_protein.hpp"
#include "zeroptmfilter/mass_match.hpp"

namespace prot {

MassMatch::MassMatch(std::vector<std::vector<int>> &mass_2d, 
                     std::vector<std::vector<double>> &real_shift_2d,
                     std::vector<std::vector<int>> &pos_2d,
                     double max_proteoform_mass, double scale) {
  scale_ = scale;
  LOG_DEBUG("Scale: " << scale_);
  LOG_DEBUG("Proteoform number: " << mass_2d.size());

  col_num_ = max_proteoform_mass * scale_;
  proteo_num_ = mass_2d.size();
  LOG_DEBUG("column number: " << col_num_);

  LOG_DEBUG("start init");
  initProteoformBeginEnds(mass_2d, real_shift_2d);
  LOG_DEBUG("row number: " << row_num_);

  LOG_DEBUG("init indexes");
  initIndexes(mass_2d, real_shift_2d, pos_2d);
}

MassMatch::~MassMatch(){
}

inline void MassMatch::initProteoformBeginEnds(std::vector<std::vector<int>> &mass_2d,
                                               std::vector<std::vector<double>> &shift_2d) {
  //no need to init
  proteo_row_begins_.resize(proteo_num_);
  proteo_row_ends_.resize(proteo_num_);
  int  pnt = 0;
  for(int i=0; i< proteo_num_; i++){
    proteo_row_begins_[i] = pnt;
    int len = shift_2d[i].size();
    proteo_row_ends_[i] = pnt + len - 1;
    //LOG_DEBUG("begin " << proteo_row_begins_[i] << " end " << proteo_row_ends_[i]);
    pnt += len;
  }
  row_num_ = pnt;
  row_proteo_ids_.resize(row_num_);
  trunc_shifts_.resize(row_num_);
  for(int i =0; i < proteo_num_; i++){
    for(int j= proteo_row_begins_[i]; j<= proteo_row_ends_[i];j++){
      row_proteo_ids_[j] = i;
      int pos = j - proteo_row_ends_[i];
      trunc_shifts_[j] = shift_2d[i][pos];
    }
  }
}

inline std::vector<std::vector<int>> convertToInt(
    std::vector<std::vector<double>> &mass_2d, double scale) {
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

inline void MassMatch::compColumnMatchNums(std::vector<std::vector<int>> &mass_2d,
                                           std::vector<std::vector<int>> &shift_2d,     
                                           std::vector<std::vector<int>> &pos_2d,     
                                           std::vector<int> &col_match_nums) {
  for (size_t i = 0; i < mass_2d.size(); i++) {
    for (size_t s = 0; s < shift_2d[i].size(); s++)  {
      for (size_t cur = pos_2d[i][s]; cur < mass_2d[i].size(); cur++) {
        int shift_mass = mass_2d[i][cur] + shift_2d[i][s];
        if (shift_mass > 0) {
          if (shift_mass < col_num_) {
            /*
            if (shift_mass >= 545694 && shift_mass <= 545712) {
              LOG_DEBUG("i " << i << " s " << s << " cur " << cur << " cur mass " << mass_2d[i][cur] << " shift mass " << shift_2d[i][s]);
            }
            */
            col_match_nums[shift_mass]++;
          }
          else {
            break;
          }
        }
      }
    }
  }
}

inline void MassMatch::fillColumnIndex(std::vector<std::vector<int>> &mass_2d,
                                       std::vector<std::vector<int>> &shift_2d,     
                                       std::vector<std::vector<int>> &pos_2d,     
                                       std::vector<int> &col_index_pnts) {
  for (size_t i = 0; i < mass_2d.size(); i++) {
    if (i/1000*1000 == i) {
      LOG_DEBUG("preprocessing proteoform " << i);
    }
    for (size_t s = 0; s < shift_2d[i].size(); s++)  {
      for (size_t cur = pos_2d[i][s]; cur < mass_2d[i].size(); cur++) {
        int shift_mass = mass_2d[i][cur] + shift_2d[i][s];
        /*
        if (s == 0 || s == 1) {
          LOG_DEBUG("i " << i << " s " << s << " cur " << cur << 
                    " shift_mass " << shift_mass << " col num " << col_num_);
        }
        */
        if (shift_mass > 0) {
          if (shift_mass < col_num_) {
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
}


inline void MassMatch::initIndexes(std::vector<std::vector<int>> &mass_2d,
                                   std::vector<std::vector<double>> &real_shift_2d,
                                   std::vector<std::vector<int>> &pos_2d) {     

  std::vector<std::vector<int>> shift_2d = convertToInt(real_shift_2d, scale_);
  LOG_DEBUG("column num " << col_num_);
  std::vector<int> col_match_nums (col_num_, 0); 
  // no need to initalize 
  std::vector<int> col_index_pnts (col_num_);
  col_index_begins_.resize(col_num_);
  col_index_ends_.resize(col_num_);

  compColumnMatchNums(mass_2d, shift_2d, pos_2d, col_match_nums);

  int pnt = 0;
  for(int i=0; i< col_num_; i++){
    col_index_begins_[i] = pnt;
    col_index_pnts[i] = pnt;
    col_index_ends_[i] = pnt + col_match_nums[i]-1;
    /*
    if (i>= 1385237 && i <= 1385279) {
      LOG_DEBUG(i << " bgn " << col_index_begins_[i] << " end " << col_index_ends_[i]);
    }
    */
    pnt += col_match_nums[i];
  }
  // no need to init
  col_indexes_.resize(pnt, 0);
  LOG_DEBUG("indexes size: "<< pnt);
  fillColumnIndex(mass_2d, shift_2d, pos_2d, col_index_pnts);

}

void MassMatch::compScores(const std::vector<std::pair<int,int>> &pref_mass_errors,
                           std::vector<short> &scores) {
  compScores(pref_mass_errors, 0, 0.0, scores);
}

void MassMatch::compScores(const std::vector<std::pair<int,int>> &pref_mass_errors,
                           int start, double shift, std::vector<short> &scores) {
  int begin_index;
  int end_index;
  int m;
  for(size_t i = start; i<pref_mass_errors.size(); i++){
    m = pref_mass_errors[i].first + shift;
    // m - errors[i] performs better than m - errors[i] -  errors[bgn_pos]
    int left = m - pref_mass_errors[i].second;
    if(left < 0){
      left=0;
    }
    int right = m + pref_mass_errors[i].second;
    if(right < 0 || right >= col_num_){
      continue;
    }
    //LOG_DEBUG("left " << left << " right " << right);
    begin_index = col_index_begins_[left];
    end_index= col_index_ends_[right];
    for(int j=begin_index;j<=end_index;j++){
      scores[col_indexes_[j]]++;
      //LOG_DEBUG("ROW INDEX " << col_indexes_[j] << " m " << m << " left " << left << " right " << right << " score " << scores[col_indexes_[j]]);
    }
  }
}

void MassMatch::compMatchScores(const std::vector<std::pair<int,int>> &pref_mass_errors,
                                const std::pair<int,int> &prec_minus_water_mass_error, 
                                std::vector<short> &scores) {
  compScores(pref_mass_errors, 0, 0.0, scores);
  // precursor mass 
  int begin_index, end_index;
  int m = prec_minus_water_mass_error.first;
  // m - errors[i] performs better than m - errors[i] -  errors[bgn_pos]
  int left = m - prec_minus_water_mass_error.second;
  if(left < 0){
    left=0;
  }
  int right = m + prec_minus_water_mass_error.second;
  if(right >= 0 && right < col_num_){
    // update scores
    begin_index = col_index_begins_[left];
    end_index= col_index_ends_[right];
    //LOG_DEBUG("prec left " << left << " pref right " << right 
    //          << " begin index " << begin_index << " end index " << end_index);
    for(int j=begin_index;j<=end_index;j++){
      scores[col_indexes_[j]] += getPrecursorMatchScore();
    }
  }
}

} /* namespace prot */
