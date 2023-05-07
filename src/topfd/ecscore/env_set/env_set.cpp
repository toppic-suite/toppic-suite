//Copyright (c) 2014 - 2022, The Trustees of Indiana University.
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

#include "topfd/ecscore/env_set/env_set.hpp"

namespace toppic {

EnvSet::EnvSet(const SeedEnvelopePtr seed_ptr, ExpEnvelopePtrVec env_list, 
               int start, int end, double noise_inte, double sn_ratio) {
  seed_ptr_ = seed_ptr;
  exp_env_list_ = env_list;
  start_spec_id_ = start;
  end_spec_id_ = end;
  initMedianXic(noise_inte, sn_ratio);
}

void EnvSet::initMedianXic(double noise_inte, double snr) {
  std::vector<double> theo_envelope_inte = seed_ptr_->getInteList();
  int refer_idx = std::max_element(theo_envelope_inte.begin(), theo_envelope_inte.end()) 
                  - theo_envelope_inte.begin();
  std::vector<double> inte_list;
  std::vector<double> env_inte_list;
  for (ExpEnvelopePtr env: exp_env_list_) {
    if (env->getMatchPeakNum(refer_idx) < 2) {
      inte_list.push_back(0);
      env_inte_list.push_back(0);
      continue;
    }
    std::vector<double> exp_envelope_inte = env->getInteList();
    /*
    double ratio = env_utils::calcInteRatio_scan(theo_envelope_inte, exp_envelope_inte);
    double theoretical_peak_sum = 0;
    theoretical_peak_sum = theoretical_peak_sum + (theo_envelope_inte[refer_idx] * ratio);
    if (refer_idx - 1 >= 0)
      theoretical_peak_sum = theoretical_peak_sum + (theo_envelope_inte[refer_idx - 1] * ratio);
    if (refer_idx + 1 < static_cast<int>(theo_envelope_inte.size()))
      theoretical_peak_sum = theoretical_peak_sum + (theo_envelope_inte[refer_idx + 1] * ratio);
    inte_list.push_back(theoretical_peak_sum);

    std::vector<double> scaled_theo_inte;
    for (auto inte: theo_envelope_inte) {
      double scaled_inte = inte * ratio;
      if (scaled_inte >= (noise_inte_level * snr))
        scaled_theo_inte.push_back(scaled_inte);
      else
        scaled_theo_inte.push_back(0);
    }
    env_inte_list.push_back(std::accumulate(scaled_theo_inte.begin(), scaled_theo_inte.end(), 0));
    */
  }
  /*
  Xic xic = Xic(start_spec_id_, seed_env_.getSpecId(), inte_list, env_inte_list);
  */
}

void EnvSet::setSpecId(int start_spec_id, int end_spec_id) {
  start_spec_id_ = start_spec_id;
  end_spec_id_ = end_spec_id;
}

}

/**

void EnvSet::get_weight_mz_error(double *weight_sum, double *error_sum) {
  std::vector<SimplePeak> seed_peak_list = seed_env_.getPeakList();
  std::vector<double> inte_list = xic_.getInteList();
  int num_exp_env = exp_env_list_.size();
  int num_peaks_theo_env = seed_peak_list.size();
  *weight_sum = 0;
  *error_sum = 0;
  for (int exp_env_id = 0; exp_env_id < num_exp_env; exp_env_id++) {
    ExpEnvelope expEnvelope = exp_env_list_[exp_env_id];
    if (expEnvelope.isEmpty())
      continue;
    std::vector<ExpPeak> exp_peak_list = expEnvelope.getExpEnvList();
    for (int peak_idx = 0; peak_idx < num_peaks_theo_env; peak_idx++) {
      ExpPeak peak = exp_peak_list[peak_idx];
      if (!peak.isEmpty()) {
        double cur_inte = seed_peak_list[peak_idx].getInte() * inte_list[exp_env_id];
        double cur_err = peak.getPos() - seed_peak_list[peak_idx].getPos();
        *error_sum = *error_sum + (cur_inte * cur_err);
        *weight_sum = *weight_sum + cur_inte;
      }
    }
  }
}

std::vector<double> EnvSet::comp_exp_inte_sum_list() {
  int peak_num = seed_env_.get_peak_num();
  std::vector<double> sum_list(peak_num, 0);
  for (auto &env: exp_env_list_) {
    if (env.isEmpty())
      continue;
    std::vector<double> cur_sum_list = env.get_inte_list();
    for (int p_i = 0; p_i < peak_num; p_i++)
      sum_list[p_i] = sum_list[p_i] + cur_sum_list[p_i];
  }
  return sum_list;
}

void EnvSet::shortlistExpEnvs() {
  std::vector<double> inte_list = xic_.getInteList();
  std::vector<double> smoothed_inte_list = xic_.getSmoothedInteList();
  std::vector<double> env_inte_list = xic_.getEnvInteList();

  std::vector<double> shortlisted_inte_list;
  std::vector<double> shortlisted_smoothed_inte_list;
  std::vector<double> shortlisted_env_inte_list;
  int num_exp_env = exp_env_list_.size();
  std::vector<ExpEnvelope> tmp;
  for (int i = 0; i < num_exp_env; i++) {
    if (exp_env_list_[i].getSpecId() >= start_spec_id_ and exp_env_list_[i].getSpecId() <= end_spec_id_) {
      tmp.push_back(exp_env_list_[i]);
      shortlisted_inte_list.push_back(inte_list[i]);
      shortlisted_smoothed_inte_list.push_back(smoothed_inte_list[i]);
      shortlisted_env_inte_list.push_back(env_inte_list[i]);
    }
  }
  exp_env_list_ = tmp;
  xic_ = Xic(exp_env_list_[0].getSpecId(), seed_env_.getSpecId(), shortlisted_inte_list,
             shortlisted_smoothed_inte_list, shortlisted_env_inte_list);
}

double _get_left_max(int pos, std::vector<double> &y) {
  double max_val = -100000000;
  for (int i = 0; i < pos; i++) {
    if (y[i] > max_val)
      max_val = y[i];
  }
  return max_val;
}

double _get_right_max(int pos, std::vector<double> &y) {
  double max_val = -100000000;
  int vec_length = y.size();
  for (int i = pos + 1; i < vec_length; i++) {
    if (y[i] > max_val)
      max_val = y[i];
  }
  return max_val;
}

void EnvSet::refine_feature_boundary() {
  double split_feature_intensity_ratio = 0.4;
  int base_spec = this->seed_env_.getSpecId() - this->start_spec_id_;
  std::vector<double> env_xic = this->xic_.getInteList();
  std::vector<double> smoothed_env_xic = this->xic_.getSmoothedInteList();

  /// Left side
  std::vector<double> left_data(smoothed_env_xic.begin(), smoothed_env_xic.begin() + base_spec + 1);
  std::vector<int> minima_left = utility_functions::findLocalMinima(left_data);
  std::vector<double> minima_vals_left;
  for (auto m: minima_left) minima_vals_left.push_back(left_data[m]);
  int start_split_point = this->start_spec_id_;
  while (!minima_vals_left.empty()) {
    int idx = std::min_element(minima_vals_left.begin(), minima_vals_left.end()) - minima_vals_left.begin();
    int pos = minima_left[idx];
    minima_vals_left.erase(minima_vals_left.begin() + idx);
    double leftMax = _get_left_max(pos, left_data);
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
  std::vector<double> right_data(smoothed_env_xic.begin() + base_spec, smoothed_env_xic.end());
  std::vector<int> minima_right = utility_functions::findLocalMinima(right_data);
  std::vector<double> minima_vals_right;
  for (auto m: minima_right) minima_vals_right.push_back(right_data[m]);
  int end_split_point = -1;
  while (!minima_vals_right.empty()) {
    int idx = std::min_element(minima_vals_right.begin(), minima_vals_right.end()) - minima_vals_right.begin();
    int pos = minima_right[idx];
    minima_vals_right.erase(minima_vals_right.begin() + idx);
    double rightMax = _get_right_max(pos, right_data);
    if (rightMax == 0) continue;
    if (right_data[pos] / rightMax <= split_feature_intensity_ratio) {
      end_split_point = pos;
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

  this->setSpecId(start, end);
  this->shortlistExpEnvs();
}

std::vector<std::vector<double>> EnvSet::get_map(double snr, double noise_inte) {
  std::vector<std::vector<double>> map;
  for (auto &exp_env: exp_env_list_) {
    std::vector<double> exp_peaks = exp_env.get_inte_list();
    std::vector<double> theo_peaks = seed_env_.get_inte_list();
    double inte_ratio = env_utils::calcInteRatio_scan(theo_peaks, exp_peaks);
    std::vector<double> scalled_theo_env;
    for (double theo_peak: theo_peaks) {
      double scaled_inte = theo_peak * inte_ratio;
      if (scaled_inte > snr * noise_inte)
        scalled_theo_env.push_back(scaled_inte);
      else
        scalled_theo_env.push_back(0);
    }
    map.push_back(scalled_theo_env);
  }
  return map;
}

double EnvSet::comp_intensity(double snr, double noise_inte) {
  std::vector<std::vector<double>> map = this->get_map(snr, noise_inte);
  if (map.size() == 0)
    return 0;
  std::vector<double> aggregate_inte(map[0].size(), 0.0);
  int num_peaks = aggregate_inte.size();
  for (auto &sp_map: map)
    for (int peakIdx = 0; peakIdx < num_peaks; peakIdx++)
      aggregate_inte[peakIdx] = aggregate_inte[peakIdx] + sp_map[peakIdx];
  double abundance = std::accumulate(aggregate_inte.begin(), aggregate_inte.end(), 0.0);
  return abundance;
}

void EnvSet::remove_peak_data(PeakMatrix &peakMatrix) {
  std::vector<std::vector<double>> map = get_map(3.0, peakMatrix.get_min_inte());
  int num_exp_env = exp_env_list_.size();
  for (int env_id = 0; env_id < num_exp_env; env_id++) {
    ExpEnvelope exp_env = exp_env_list_[env_id];
    int spec_id = exp_env.getSpecId();
    if (spec_id < 0 or spec_id >= peakMatrix.get_spec_num())
      continue;
    std::vector<ExpPeak> exp_data = exp_env.getExpEnvList();
    std::vector<double> theo_data = map[env_id];
    int num_peaks_exp_env = exp_data.size();
    for (int peak_id = 0; peak_id < num_peaks_exp_env; peak_id++) {
      ExpPeak exp_peak = exp_data[peak_id];
      if (exp_peak.isEmpty())
        continue;
      int bin_idx = peakMatrix.get_index(exp_peak.getPos());
      std::vector<ExpPeak> peaks = peakMatrix.get_bin_peak(spec_id, bin_idx);
      double theo_peak = theo_data[peak_id];
      std::vector<ExpPeak> other_peaks;
      for (auto peak: peaks) {
        if (std::abs(peak.getInte() - exp_peak.getInte()) < 0.0001) {
          if (peak.getInte() / theo_peak < 4)
            continue;
          else {
            if (peak.getInte() - theo_peak > 0)
              peak.setInte(peak.getInte() - theo_peak);
            else
              peak.setInte(0);
          }
        }
        other_peaks.push_back(peak);
      }
      peakMatrix.set_bin_peak(spec_id, bin_idx, other_peaks);
    }
  }
}
**/


