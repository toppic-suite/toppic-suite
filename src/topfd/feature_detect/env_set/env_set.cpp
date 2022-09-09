//
// Created by abbash on 7/29/22.
//

#include "env_set.hpp"
#include "topfd/feature_detect/score/env_util.hpp"
#include "topfd/feature_detect/util/utility_functions.hpp"
#include <numeric>
#include <algorithm>
#include <iostream>

toppic::EnvSet::EnvSet() {
  seed_env_ = SeedEnvelope();
  start_spec_id_ = -1;
  end_spec_id_ = -1;
  xic_ = Xic();
}

toppic::EnvSet::EnvSet(const EnvSet & es){
  seed_env_ = es.seed_env_;
  for (const auto& exp_env : es.exp_env_list_)
    exp_env_list_.push_back(exp_env);
  start_spec_id_ = es.start_spec_id_;
  end_spec_id_ = es.end_spec_id_;
  xic_ = Xic(es.xic_);
}

toppic::EnvSet::EnvSet(const SeedEnvelope& envelope, std::vector<ExpEnvelope> env_list, int start, int end) {
  seed_env_ = envelope;
  exp_env_list_ = env_list;
  start_spec_id_ = start;
  end_spec_id_ = end;
  xic_ = init_median_xic();
}

//toppic::Xic toppic::EnvSet::init_median_xic(){
//  std::vector<double> inte_list;
//  for (const ExpEnvelope& env : exp_env_list_){
//    double ratio = get_median_ratio(env);
//    inte_list.push_back(ratio);
//  }
//  Xic xic = Xic(start_spec_id_, seed_env_.getSpecId(), inte_list);
//  return xic;
//}

toppic::Xic toppic::EnvSet::init_median_xic(){
  std::vector<double> theo_envelope_inte = seed_env_.get_inte_list();
  int refer_idx = std::max_element(theo_envelope_inte.begin(), theo_envelope_inte.end()) - theo_envelope_inte.begin();

  std::vector<double> inte_list;
  for (ExpEnvelope& env : exp_env_list_){
    if (env.get_match_peak_num(refer_idx) < 2) {
      inte_list.push_back(0);
      continue;
    }
    double ratio = env_utils::calcInteRatio_scan(theo_envelope_inte, env.get_inte_list());
    double theoretical_peak_sum = 0;
    theoretical_peak_sum = theoretical_peak_sum + (theo_envelope_inte[refer_idx] * ratio);
    if (refer_idx - 1 >= 0)
      theoretical_peak_sum = theoretical_peak_sum + (theo_envelope_inte[refer_idx-1] * ratio);
    if (refer_idx + 1 < theo_envelope_inte.size())
      theoretical_peak_sum = theoretical_peak_sum + (theo_envelope_inte[refer_idx+1] * ratio);
    inte_list.push_back(theoretical_peak_sum);
  }
  Xic xic = Xic(start_spec_id_, seed_env_.getSpecId(), inte_list);
  return xic;
}

double find_median(std::vector<double> in_data) {
  int n = in_data.size();
  if (n % 2 == 0) {
    std::nth_element(in_data.begin(), in_data.begin() + n / 2, in_data.end());
    nth_element(in_data.begin(), in_data.begin() + (n - 1) / 2, in_data.end());
    return (double)(in_data[(n - 1) / 2] + in_data[n / 2]) / 2.0;
  }
  else {
    nth_element(in_data.begin(), in_data.begin() + n / 2, in_data.end());
    return (double)in_data[n / 2];
  }
}

double toppic::EnvSet::get_median_ratio(ExpEnvelope env){
  return env_utils::calcInteRatio_scan(seed_env_.get_inte_list(), env.get_inte_list());
}


//double toppic::EnvSet::get_median_ratio(ExpEnvelope env){
//  std::vector<ExpPeak> exp_peak_list = env.getExpEnvList();
//  std::vector<SimplePeak> theo_peak_list = seed_env_.getPeakList();
//  std::vector<double> ratio_list;
//  for (int i = 0; i < theo_peak_list.size(); i++){
//    double ratio;
//    if (!exp_peak_list[i].isEmpty()) {
//      ratio = exp_peak_list[i].getInte() / theo_peak_list[i].getInte();
//      ratio_list.push_back(ratio);
//    }
//  }
//  double ratio_sum = std::accumulate(ratio_list.begin(), ratio_list.end(), 0.0);
//  if (ratio_sum == 0)
//    return 0;
//  else
//    return find_median(ratio_list);
//}

void toppic::EnvSet::get_coordinates(spec_list spectra_list, std::vector<double> x, std::vector<double> y, std::vector<double> z){
  for (auto env : exp_env_list_) {
    for (int peak_id = 0; peak_id < env.get_peak_num(); peak_id++){
      ExpPeak peak = env.get_peak(peak_id);
      if (!peak.isEmpty()) {
        int spec_id = env.getSpecId();
        x.push_back(peak.getPos());
        y.push_back(spectra_list[spec_id].getRt());
        z.push_back(peak.getInte());
      }
    }
  }
}

double toppic::EnvSet::get_seed_inte_ratio(){
  int idx = seed_env_.getSpecId() - start_spec_id_;
  double ratio = get_median_ratio(exp_env_list_[idx]);
  return ratio;
}

void toppic::EnvSet::remove_non_consecutive_peaks(int i, int max_miss_peak){
  //search backward
  int idx = seed_env_.getSpecId() - start_spec_id_;
  int miss_num = 0;
  while (idx >= 0){
    std::vector<ExpPeak> peaks = exp_env_list_[idx].getExpEnvList();
    if (peaks[i].isEmpty())
      miss_num = miss_num + 1;
    else
      miss_num = 0;
    if (miss_num >= max_miss_peak)
      break;
    idx = idx - 1;
  }
  while (idx >= 0) {
    std::vector<ExpPeak> peaks = exp_env_list_[idx].getExpEnvList();
    peaks[i] = ExpPeak();
    exp_env_list_[idx].setExpEnvList(peaks);
    idx = idx - 1;
  }

  // search forward
  idx = seed_env_.getSpecId() - start_spec_id_;
  miss_num = 0;
  while (idx < exp_env_list_.size()) {
    std::vector<ExpPeak> peaks = exp_env_list_[idx].getExpEnvList();
    if (peaks[i].isEmpty())
      miss_num = miss_num + 1;
    else
      miss_num = 0;
    if (miss_num >= max_miss_peak)
      break;
    idx = idx + 1;
  }
  while (idx < exp_env_list_.size()) {
    std::vector<ExpPeak> peaks = exp_env_list_[idx].getExpEnvList();
    peaks[i] = ExpPeak();
    exp_env_list_[idx].setExpEnvList(peaks);
    idx = idx + 1;
  }
}

void toppic::EnvSet::remove_all_non_consecutive_peaks(int max_miss_peak){
  for (int i = 0; i < seed_env_.get_peak_num(); i++)
    remove_non_consecutive_peaks(i, max_miss_peak);
}

void toppic::EnvSet::simple_remove_matrix_peaks(PeakMatrix peak_matrix){
  for (auto env : exp_env_list_) {
    for (int peak_id = 0; peak_id < env.get_peak_num(); peak_id++) {
      ExpPeak peak = env.get_peak(peak_id);
      if (!peak.isEmpty())
        peak_matrix.remove_peak(peak);
    }
  }
}

void toppic::EnvSet::remove_matrix_peaks(PeakMatrix& peak_matrix){
  for (auto env : exp_env_list_) {
    int peak_num = 0;
    for (int peak_id = 0; peak_id < env.get_peak_num(); peak_id++) {
      ExpPeak peak = env.get_peak(peak_id);
      if (!peak.isEmpty()){
        peak_num = peak_num + 1;
        peak_matrix.remove_peak(peak);
      }
    }
    if (peak_num > 0) {
      int spec_id = env.getSpecId();
      double min_pos = 0, max_pos = 0;
      env.get_min_max_pos(&min_pos, &max_pos);
      peak_matrix.remove_peak_in_range(spec_id, min_pos, max_pos);
    }
  }
}

void toppic::EnvSet::get_weight_mz_error(double* weight_sum, double* error_sum){
  std::vector<SimplePeak> seed_peak_list = seed_env_.getPeakList();
  *weight_sum = 0;
  *error_sum = 0;
  for (int spec_idx = 0; spec_idx < exp_env_list_.size(); spec_idx++) {
    ExpEnvelope expEnvelope = exp_env_list_[spec_idx];
    if (expEnvelope.isEmpty())
      continue;
    std::vector<ExpPeak> exp_peak_list = expEnvelope.getExpEnvList();
    for (int peak_idx = 0; peak_idx < seed_peak_list.size(); peak_idx++){
      ExpPeak peak = exp_peak_list[peak_idx];
      if (!peak.isEmpty()) {
        std::vector<double> inte_list = xic_.getInteList();
        double cur_inte = seed_peak_list[peak_idx].getInte() * inte_list[spec_idx];
        double cur_err = peak.getPos() - seed_peak_list[peak_idx].getPos();
        *error_sum = *error_sum + (cur_inte * cur_err);
        *weight_sum = *weight_sum + cur_inte;
      }
    }
  }
}

std::vector<double> toppic::EnvSet::comp_exp_inte_sum_list() {
  int peak_num = seed_env_.get_peak_num();
  std::vector<double> sum_list (peak_num, 0);
  for (auto env: exp_env_list_){
    if (env.isEmpty())
      continue;
    std::vector<double> cur_sum_list = env.get_inte_list();
    for (int p_i = 0; p_i < peak_num; p_i++)
      sum_list[p_i] = sum_list[p_i] + cur_sum_list[p_i];
  }
  return sum_list;
}

double get_left_max(int pos, std::vector<double> y) {
  double max_val = -100000000;
  for (int i = 0; i < pos; i++) {
    if (y[i] > max_val)
      max_val = y[i];
  }
  return max_val;
}

double get_right_max(int pos, std::vector<double> y){
  double max_val = -100000000;
  for (int i = pos + 1; i < y.size(); i++) {
    if (y[i] > max_val)
      max_val = y[i];
  }
  return max_val;
}

void toppic::EnvSet::refine_feature_boundary(){
  double split_feature_intensity_ratio = 0.4;
  int base_spec = this->seed_env_.getSpecId() - this->start_spec_id_;
  std::vector<double> env_xic = this->xic_.getInteList();
  std::vector<double> smoothed_env_xic = this->xic_.getSmoothedInteList();
//  for (int i = 0; i < env_xic.size(); i++)    std::cout << i << ", " << env_xic[i] << ", " << smoothed_env_xic[i] << std::endl;

  /// Left side
  std::vector<double> left_data(smoothed_env_xic.begin(), smoothed_env_xic.begin()+base_spec+1);
//  for (auto l : left_data) std::cout << "Left: " << l << std::endl;
  std::vector<double> minima_left = utility_functions::findLocalMinima(left_data);
  std::vector<double> minima_vals_left;
  for (auto m: minima_left) minima_vals_left.push_back(left_data[m]);
//  for (int i = 0; i < minima_left.size(); i++) std::cout << "left_min: " << minima_left[i] << ", " << minima_vals_left[i] << std::endl;
  int start_split_point = this->start_spec_id_;
  while (minima_vals_left.size() > 0) {
    int idx = std::min_element(minima_vals_left.begin(), minima_vals_left.end()) - minima_vals_left.begin();
    int pos = minima_left[idx];
    minima_vals_left.erase(minima_vals_left.begin() + idx);
    double leftMax = get_left_max(pos, left_data);
    if (leftMax == 0) continue;
    if (left_data[pos] / leftMax <= split_feature_intensity_ratio) {
      start_split_point = start_split_point + pos;
      std::vector<double> temp_left_data(left_data.begin() + pos, left_data.end());
      left_data = temp_left_data;
      minima_left = utility_functions::findLocalMinima(left_data);
      minima_vals_left.clear();
      for (auto m: minima_left) minima_vals_left.push_back(left_data[m]);
    }
  }

  /// Right side
  std::vector<double> right_data(smoothed_env_xic.begin()+base_spec, smoothed_env_xic.end());
//  for (auto r : right_data) std::cout << "Right: " << r << std::endl;
  std::vector<double> minima_right = utility_functions::findLocalMinima(right_data);
  std::vector<double> minima_vals_right;
  for (auto m: minima_right) minima_vals_right.push_back(right_data[m]);
//  for (int i = 0; i < minima_right.size(); i++) std::cout << "Right_min: " << minima_right[i] << ", " << minima_vals_right[i] << std::endl;
  int end_split_point = -1;
  while (minima_vals_right.size() > 0) {
    int idx = std::min_element(minima_vals_right.begin(), minima_vals_right.end()) - minima_vals_right.begin();
    int pos = minima_right[idx];
    minima_vals_right.erase(minima_vals_right.begin() + idx);
    double rightMax = get_right_max(pos, right_data);
    if (rightMax == 0) continue;
//    std::cout << idx << ", " << pos << ", " << rightMax << ", " << right_data[pos] << ", " << right_data[pos]/rightMax << std::endl;
    if (right_data[pos] / rightMax <= split_feature_intensity_ratio) {
      end_split_point =  pos;
//      std::cout << "SPLITTTTT: " <<  end_split_point << std::endl;
      std::vector<double> temp_right_data(right_data.begin(), right_data.begin() + pos - 1);
      right_data = temp_right_data;
      minima_right = utility_functions::findLocalMinima(right_data);
      minima_vals_right.clear();
      for (auto m: minima_right) minima_vals_right.push_back(right_data[m]);
    }
  }
  int start = this->start_spec_id_;
  if (start_split_point > -1)
    start = start_split_point;
  int end = this->end_spec_id_;
  if (end_split_point > -1)
    end = seed_env_.getSpecId() + end_split_point;

//  std::cout << "New Boundary: " << start << ", " << end << std::endl;
  this->setStartSpecId(start);
  this->setEndSpecId(end);
  this->xic_.setStartSpecId(start);
}

std::vector<std::vector<double>> toppic::EnvSet::get_map(double snr, double noise_inte){
  std::vector<std::vector<double>> map;
  std::vector<ExpEnvelope> exp_envs = this->exp_env_list_;
  for (auto exp_env: exp_envs) {
    std::vector<double> exp_peaks = exp_env.get_inte_list();
    std::vector<double> theo_peaks = this->get_theo_distribution_inte();
    double inte_ratio = env_utils::calcInteRatio_scan(theo_peaks, exp_peaks);
    std::vector<double> scalled_theo_env;
    for (int i = 0; i < theo_peaks.size(); i++) {
      double scaled_inte = theo_peaks[i] * inte_ratio;
      if (scaled_inte > snr * noise_inte)
        scalled_theo_env.push_back(scaled_inte);
      else
        scalled_theo_env.push_back(0);
    }
    map.push_back(scalled_theo_env);
  }
  return map;
}

double toppic::EnvSet::comp_intensity(double snr, double noise_inte){
  std::vector<std::vector<double>> map = this->get_map(snr, noise_inte);
  std::vector<double> aggregate_inte (map[0].size(), 0.0);
  for (int spId = 0; spId < map.size(); spId++)
    for (int peakIdx = 0; peakIdx < aggregate_inte.size(); peakIdx++)
        aggregate_inte[peakIdx] = aggregate_inte[peakIdx] + map[spId][peakIdx];
//  for (auto i : aggregate_inte) std::cout << "Abundance - aggregate inte: " << i << std::endl;
  double abundance = std::accumulate(aggregate_inte.begin(), aggregate_inte.end(), 0.0);
//  std::cout << "Abundance: " << abundance << std::endl;
  return abundance;
}