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

#include <numeric>

#include "common/util/logger.hpp"

#include "topfd/ecscore/envelope/env_util.hpp"
#include "topfd/ecscore/env_set/env_set_util.hpp"
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

void EnvSet::initMedianXic(double noise_inte, double sn_ratio) {
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
    double ratio = env_set_util::calcInteRatio(theo_envelope_inte, exp_envelope_inte);
    double top_three_inte = theo_envelope_inte[refer_idx] * ratio;
    if (refer_idx - 1 >= 0)
      top_three_inte += theo_envelope_inte[refer_idx - 1] * ratio;
    if (refer_idx + 1 < static_cast<int>(theo_envelope_inte.size()))
      top_three_inte += theo_envelope_inte[refer_idx + 1] * ratio;
    inte_list.push_back(top_three_inte);

    std::vector<double> scaled_theo_inte;
    for (auto inte: theo_envelope_inte) {
      double scaled_inte = inte * ratio;
      if (scaled_inte >= (noise_inte * sn_ratio))
        scaled_theo_inte.push_back(scaled_inte);
      else
        scaled_theo_inte.push_back(0);
    }
    env_inte_list.push_back(std::accumulate(scaled_theo_inte.begin(), scaled_theo_inte.end(), 0));
  }
  xic_ptr_ = std::make_shared<Xic>(start_spec_id_, seed_ptr_->getSpecId(), inte_list, env_inte_list);
}

void EnvSet::setSpecId(int start_spec_id, int end_spec_id) {
  start_spec_id_ = start_spec_id;
  end_spec_id_ = end_spec_id;
}

std::vector<double> EnvSet::compExpInteSumList() {
  int peak_num = seed_ptr_->getPeakNum();
  std::vector<double> sum_list(peak_num, 0);
  for (auto &env: exp_env_list_) {
    if (env == nullptr) {
      continue;
    }
    std::vector<double> cur_sum_list = env->getInteList();
    for (int p_i = 0; p_i < peak_num; p_i++)
      sum_list[p_i] = sum_list[p_i] + cur_sum_list[p_i];
  }
  return sum_list;
}

void EnvSet::getWeightMzError(double &weight_sum, double &error_sum) {
  EnvSimplePeakPtrVec seed_peak_list = seed_ptr_->getPeakList();
  std::vector<double> inte_list = xic_ptr_->getInteList();
  int num_exp_env = exp_env_list_.size();
  int num_peaks_theo_env = seed_peak_list.size();
  weight_sum = 0;
  error_sum = 0;
  for (int exp_env_id = 0; exp_env_id < num_exp_env; exp_env_id++) {
    ExpEnvelopePtr env_ptr = exp_env_list_[exp_env_id];
    if (env_ptr == nullptr) {
      continue;
    }
    MatrixPeakPtrVec exp_peak_list = env_ptr->getExpPeakList();
    for (int peak_idx = 0; peak_idx < num_peaks_theo_env; peak_idx++) {
      MatrixPeakPtr peak = exp_peak_list[peak_idx];
      if (peak != nullptr) {
        double cur_inte = seed_peak_list[peak_idx]->getIntensity() * inte_list[exp_env_id];
        double cur_err = peak->getPosition() - seed_peak_list[peak_idx]->getPosition();
        error_sum = error_sum + (cur_inte * cur_err);
        weight_sum = weight_sum + cur_inte;
      }
    }
  }
}

std::vector<std::vector<double>> EnvSet::getScaledTheoIntes(double sn_ratio, 
                                                            double noise_inte) {
  std::vector<std::vector<double>> intes;
  for (auto exp_env: exp_env_list_) {
    std::vector<double> exp_intes = exp_env->getInteList();
    std::vector<double> theo_intes = seed_ptr_->getInteList();
    double inte_ratio = env_set_util::calcInteRatio(theo_intes, exp_intes);
    std::vector<double> scalled_theo_intes;
    for (double theo_inte: theo_intes) {
      double scaled_inte = theo_inte * inte_ratio;
      if (scaled_inte > sn_ratio * noise_inte)
        scalled_theo_intes.push_back(scaled_inte);
      else
        scalled_theo_intes.push_back(0);
    }
    intes.push_back(scalled_theo_intes);
  }
  return intes;
}


double EnvSet::compIntensity(double sn_ratio, double noise_inte) {
  std::vector<std::vector<double>> intes = getScaledTheoIntes(sn_ratio, noise_inte);
  if (intes.size() == 0)
    return 0;
  int peak_num = intes[0].size();
  std::vector<double> aggregate_inte(peak_num, 0.0);
  for (auto &sp_intes: intes)
    for (int peak_idx = 0; peak_idx < peak_num; peak_idx++)
      aggregate_inte[peak_idx] = aggregate_inte[peak_idx] + sp_intes[peak_idx];
  double abundance = std::accumulate(aggregate_inte.begin(), aggregate_inte.end(), 0.0);
  return abundance;
}

double getLeftMax(int pos, std::vector<double> &y) {
  double max_val = -100000000;
  for (int i = 0; i < pos; i++) {
    if (y[i] > max_val)
      max_val = y[i];
  }
  return max_val;
}

double getRightMax(int pos, std::vector<double> &y) {
  double max_val = -100000000;
  int vec_length = y.size();
  for (int i = pos + 1; i < vec_length; i++) {
    if (y[i] > max_val)
      max_val = y[i];
  }
  return max_val;
}

void EnvSet::shortlistExpEnvs() {
  std::vector<double> inte_list = xic_ptr_->getInteList();
  std::vector<double> smoothed_inte_list = xic_ptr_->getSmoothedInteList();
  std::vector<double> env_inte_list = xic_ptr_->getEnvInteList();

  std::vector<double> shortlisted_inte_list;
  std::vector<double> shortlisted_smoothed_inte_list;
  std::vector<double> shortlisted_env_inte_list;
  int num_exp_env = exp_env_list_.size();
  ExpEnvelopePtrVec tmp;
  for (int i = 0; i < num_exp_env; i++) {
    if (exp_env_list_[i]->getSpecId() >= start_spec_id_ && 
        exp_env_list_[i]->getSpecId() <= end_spec_id_) {
      tmp.push_back(exp_env_list_[i]);
      shortlisted_inte_list.push_back(inte_list[i]);
      shortlisted_smoothed_inte_list.push_back(smoothed_inte_list[i]);
      shortlisted_env_inte_list.push_back(env_inte_list[i]);
    }
  }
  exp_env_list_ = tmp;
  xic_ptr_ = std::make_shared<Xic>(exp_env_list_[0]->getSpecId(), seed_ptr_->getSpecId(), 
                                   shortlisted_inte_list, shortlisted_smoothed_inte_list, 
                                   shortlisted_env_inte_list);
}

void EnvSet::refineFeatureBoundary() {
  double split_feature_intensity_ratio = 0.4;
  int base_spec = seed_ptr_->getSpecId() - start_spec_id_;
  std::vector<double> env_xic = xic_ptr_->getInteList();
  std::vector<double> smoothed_env_xic = xic_ptr_->getSmoothedInteList();

  /// Left side
  std::vector<double> left_data(smoothed_env_xic.begin(), smoothed_env_xic.begin() + base_spec + 1);
  std::vector<int> minima_left = env_util::findLocalMinima(left_data);
  std::vector<double> minima_vals_left;
  for (auto m: minima_left) minima_vals_left.push_back(left_data[m]);
  int start_split_point = start_spec_id_;
  while (!minima_vals_left.empty()) {
    int idx = std::min_element(minima_vals_left.begin(), minima_vals_left.end()) - minima_vals_left.begin();
    int pos = minima_left[idx];
    minima_vals_left.erase(minima_vals_left.begin() + idx);
    double left_max = getLeftMax(pos, left_data);
    if (left_max == 0) continue;
    if (left_data[pos] / left_max <= split_feature_intensity_ratio) {
      start_split_point = start_split_point + pos;
      std::vector<double> temp_left_data(left_data.begin() + pos, left_data.end());
      left_data = temp_left_data;
      minima_left = env_util::findLocalMinima(left_data);
      minima_vals_left.clear();
      for (auto m: minima_left) minima_vals_left.push_back(left_data[m]);
    }
  }

  /// Right side
  std::vector<double> right_data(smoothed_env_xic.begin() + base_spec, smoothed_env_xic.end());
  std::vector<int> minima_right = env_util::findLocalMinima(right_data);
  std::vector<double> minima_vals_right;
  for (auto m: minima_right) minima_vals_right.push_back(right_data[m]);
  int end_split_point = -1;
  while (!minima_vals_right.empty()) {
    int idx = std::min_element(minima_vals_right.begin(), minima_vals_right.end()) - minima_vals_right.begin();
    int pos = minima_right[idx];
    minima_vals_right.erase(minima_vals_right.begin() + idx);
    double right_max = getRightMax(pos, right_data);
    if (right_max == 0) continue;
    if (right_data[pos] / right_max <= split_feature_intensity_ratio) {
      end_split_point = pos;
      std::vector<double> temp_right_data(right_data.begin(), right_data.begin() + pos - 1);
      right_data = temp_right_data;
      minima_right = env_util::findLocalMinima(right_data);
      minima_vals_right.clear();
      for (auto m: minima_right) minima_vals_right.push_back(right_data[m]);
    }
  }
  int start = this->start_spec_id_;
  if (start_split_point > -1)
    start = start_split_point;
  int end = this->end_spec_id_;
  if (end_split_point > -1)
    end = seed_ptr_->getSpecId() + end_split_point;

  setSpecId(start, end);
  shortlistExpEnvs();
}

void EnvSet::removePeakData(PeakMatrixPtr matrix_ptr) {
  double sn_ratio = 3.0;
  std::vector<std::vector<double>> theo_intes = getScaledTheoIntes(sn_ratio, 
                                                                   matrix_ptr->getBaseInte());
  int num_exp_env = exp_env_list_.size();
  for (int env_id = 0; env_id < num_exp_env; env_id++) {
    LOG_DEBUG("Env id " << env_id << " " << num_exp_env);
    ExpEnvelopePtr exp_env = exp_env_list_[env_id];
    if (exp_env == nullptr) {
      continue;
    }
    int spec_id = exp_env->getSpecId();
    if (spec_id < 0 or spec_id >= matrix_ptr->getSpecNum()) {
      continue;
    }
    MatrixPeakPtrVec exp_peaks = exp_env->getExpPeakList(); 
    std::vector<double> theo_env_peak_intes = theo_intes[env_id];

    int peak_num = exp_peaks.size(); 
    for (int peak_id = 0; peak_id < peak_num; peak_id++) {
      MatrixPeakPtr exp_peak = exp_peaks[peak_id];
      LOG_DEBUG("Peak id " << peak_id << " " << peak_num);
      if (exp_peak == nullptr) {
        continue;
      }
      int bin_idx = matrix_ptr->getBinIndex(exp_peak->getPosition()); 
      MatrixPeakPtrVec bin_peaks = matrix_ptr->getBinPeakList(spec_id, bin_idx);
      double theo_peak_inte = theo_env_peak_intes[peak_id];
      MatrixPeakPtrVec remain_peaks;
      for (auto peak: bin_peaks) {
        // check if peak is the same as exp_peak
        if (peak == nullptr) {
          continue;
        }
        if (std::abs(peak->getIntensity() - exp_peak->getIntensity()) < 0.0001) {
          if (peak->getIntensity() / theo_peak_inte < 4) {
            continue;
          }
          else {
            peak->setIntensity(peak->getIntensity() - theo_peak_inte);
          }
        }
        remain_peaks.push_back(peak);
      }
      matrix_ptr->setBinPeakList(spec_id, bin_idx, remain_peaks);
    }
  }
}

}

