//
// Created by abbash on 7/29/22.
//

#include "env_set.hpp"
#include <numeric>
#include <algorithm>

toppic::EnvSet::EnvSet() {
  seed_env_ = SeedEnvelope();
  start_spec_id_ = -1;
  end_spec_id_ = -1;
  xic_ = Xic();
}

toppic::EnvSet::EnvSet(const SeedEnvelope& envelope, std::vector<ExpEnvelope> env_list, int start, int end) {
  seed_env_ = envelope;
  exp_env_list_ = env_list;
  start_spec_id_ = start;
  end_spec_id_ = end;
  xic_ = init_median_xic();
}

toppic::EnvSet::EnvSet(const EnvSet & es){
  seed_env_ = es.seed_env_;
  for (const auto& exp_env : es.exp_env_list_)
    exp_env_list_.push_back(exp_env);
  start_spec_id_ = es.start_spec_id_;
  end_spec_id_ = es.end_spec_id_;
  xic_ = Xic(es.xic_);
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

double toppic::EnvSet::get_median_ratio(ExpEnvelope env){
  std::vector<ExpPeak> exp_peak_list = env.getExpEnvList();
  std::vector<SimplePeak> theo_peak_list = seed_env_.getPeakList();
  std::vector<double> ratio_list;
  for (int i = 0; i < theo_peak_list.size(); i++){
    double ratio;
    if (!exp_peak_list[i].isEmpty()) {
      ratio = exp_peak_list[i].getInte() / theo_peak_list[i].getInte();
      ratio_list.push_back(ratio);
    }
  }
  double ratio_sum = std::accumulate(ratio_list.begin(), ratio_list.end(), 0.0);
  if (ratio_sum == 0)
    return 0;
  else
    return find_median(ratio_list);
}

toppic::Xic toppic::EnvSet::init_median_xic(){
  std::vector<double> inte_list;
  for (const ExpEnvelope& env : exp_env_list_){
    double ratio = get_median_ratio(env);
    inte_list.push_back(ratio);
  }
  Xic xic = Xic(start_spec_id_, seed_env_.getSpecId(), inte_list);
  return xic;
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

void toppic::EnvSet::remove_matrix_peaks(PeakMatrix peak_matrix){
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


