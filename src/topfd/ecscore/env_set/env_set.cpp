//Copyright (c) 2014 - 2023, The Trustees of Indiana University.
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
#include <limits>
#include <algorithm>

#include "common/util/logger.hpp"

#include "topfd/ecscore/env/ms_map_env_util.hpp"
#include "topfd/ecscore/env_set/xic_util.hpp"
#include "topfd/ecscore/env_set/env_set.hpp"

namespace toppic {

EnvSet::EnvSet(const SeedEnvPtr seed_ptr, MsMapEnvPtrVec env_list,
               int start, int end, double noise_inte, double sn_ratio) {
  seed_ptr_ = seed_ptr;
  ms_map_env_list_ = env_list;
  start_spec_id_ = start;
  end_spec_id_ = end;
  min_inte_ = noise_inte * sn_ratio;
  initMedianXic();
}

EnvSet::EnvSet(const SeedEnvPtr seed_ptr, MsMapEnvPtrVec env_list,
               int start, int end, double min_inte) {
  seed_ptr_ = seed_ptr;
  ms_map_env_list_ = env_list;
  start_spec_id_ = start;
  end_spec_id_ = end;
  min_inte_ = min_inte;
  initMedianXic();
}

void EnvSet::initMedianXic() {
  std::vector<double> inte_ratio_list;
  std::vector<double> top_three_inte_list;
  std::vector<double> all_peak_inte_list;
  for (MsMapEnvPtr env_ptr: ms_map_env_list_) {
    double ratio = ms_map_env_util::compTopThreeInteRatio(seed_ptr_, env_ptr); 
    inte_ratio_list.push_back(ratio);
    double top_three_inte = seed_ptr_->compTopThreeInteSum() * ratio;
    top_three_inte_list.push_back(top_three_inte);
    double all_peak_inte = seed_ptr_->compScaledInteSum(ratio, min_inte_);
    all_peak_inte_list.push_back(all_peak_inte); 
  }
  xic_ptr_ = std::make_shared<Xic>(inte_ratio_list, top_three_inte_list, 
                                   all_peak_inte_list);
}

double EnvSet::getXicSeedAllPeakInte() {
  int seed_spec_id = seed_ptr_->getSpecId();
  if (seed_spec_id < start_spec_id_ || seed_spec_id > end_spec_id_) {
    return 0;
  }
  int seed_idx = seed_spec_id - start_spec_id_;
  return xic_ptr_->getAllPeakInte(seed_idx);
}

int EnvSet::countEnvNum() {
  int count = 0;
  for (size_t i = 0; i < ms_map_env_list_.size(); i++) {
    if (ms_map_env_list_[i] != nullptr) {
      count++;
    }
  }
  return count;
}

std::vector<double> EnvSet::compAggrEnvInteList() {
  int peak_num = seed_ptr_->getPeakNum();
  std::vector<double> inte_list(peak_num, 0);
  for (auto &env: ms_map_env_list_) {
    if (env == nullptr) {
      continue;
    }
    std::vector<double> cur_sum_list = env->getInteList();
    for (int p_i = 0; p_i < peak_num; p_i++)
      inte_list[p_i] = inte_list[p_i] + cur_sum_list[p_i];
  }
  return inte_list;
}

std::vector<double> EnvSet::compAggrEnvMzList() {
  int peak_num = seed_ptr_->getPeakNum();
  std::vector<double> mz_list(peak_num, 0.0);
  for (int peak_idx = 0; peak_idx < peak_num; peak_idx++) {
    int weight = 0;
    for (size_t env_idx = 0; env_idx < ms_map_env_list_.size(); env_idx++) {
      MsMapPeakPtr peak_ptr = ms_map_env_list_[env_idx]->getPeakPtr(peak_idx);
      if (peak_ptr != nullptr) {
        mz_list[peak_idx] = mz_list[peak_idx]
          + (peak_ptr->getPosition() * peak_ptr->getIntensity());
        weight += peak_ptr->getIntensity();
      }
    }
    if (weight > 0) {
      mz_list[peak_idx] = mz_list[peak_idx] / weight;
    }
  }
  return mz_list; 
}


std::vector<std::vector<double>> EnvSet::getScaledTheoIntes(int min_inte) {
  std::vector<std::vector<double>> results;
  std::vector<double> ratio_list = xic_ptr_->getInteRatioList();
  for (size_t i = 0; i < ratio_list.size(); i++) {
    std::vector<double> spec_inte_list = seed_ptr_->getScaledInteList(ratio_list[i], min_inte);
    results.push_back(spec_inte_list);
  }
  return results;
}

void EnvSet::removePeakData(MsMapPtr ms_map_ptr) {
  double peak_remove_ratio = 4;
  int env_num = ms_map_env_list_.size();
  for (int env_id = 0; env_id < env_num; env_id++) {
    MsMapEnvPtr env_ptr = ms_map_env_list_[env_id];
    if (env_ptr == nullptr) {
      continue;
    }
    int spec_id = env_ptr->getSpecId();
    if (spec_id < 0 or spec_id >= ms_map_ptr->getRowNum()) {
      continue;
    }
    MsMapPeakPtrVec env_peak_list = env_ptr->getMsMapPeakList();
    double ratio = xic_ptr_->getInteRatio(env_id);
    std::vector<double> theo_env_peak_intes = seed_ptr_->getScaledInteList(ratio, min_inte_);
    int peak_num = env_peak_list.size();
    for (int peak_id = 0; peak_id < peak_num; peak_id++) {
      MsMapPeakPtr exp_peak = env_peak_list[peak_id];
      if (exp_peak == nullptr) {
        continue;
      }
      int col_idx = ms_map_ptr->getColIndex(exp_peak->getPosition());
      MsMapPeakPtrVec bin_peaks = ms_map_ptr->getBinPeakList(spec_id, col_idx);
      double theo_peak_inte = theo_env_peak_intes[peak_id];
      MsMapPeakPtrVec remain_peaks;
      for (auto peak: bin_peaks) {
        // check if peak is the same as exp_peak
        if (peak == nullptr) {
          continue;
        }
        if (peak == exp_peak) {
          // we need to double check if the parameter 4 is a good one
          if (peak->getIntensity() / theo_peak_inte < peak_remove_ratio) {
            // skip this peak
            continue;
          }
          else {
            peak->setIntensity(peak->getIntensity() - theo_peak_inte);
          }
        }
        remain_peaks.push_back(peak);
      }
      ms_map_ptr->setBinPeakList(spec_id, col_idx, remain_peaks);
    }
  }
}

void EnvSet::shortlistExpEnvs() {
  int env_num = ms_map_env_list_.size();
  MsMapEnvPtrVec tmp;
  for (int i = 0; i < env_num; i++) {
    if (ms_map_env_list_[i]->getSpecId() >= start_spec_id_ &&
        ms_map_env_list_[i]->getSpecId() <= end_spec_id_) {
      tmp.push_back(ms_map_env_list_[i]);
    }
  }
    ms_map_env_list_ = tmp;
}

std::pair<double, double> EnvSet::getMzErrorAndWeight() {
  EnvPeakPtrVec seed_peak_list = seed_ptr_->getPeakPtrList();
  std::vector<double> ratio_list = xic_ptr_->getInteRatioList();
  int num_peaks_theo_env = seed_peak_list.size();
  double weight_sum = 0;
  double error_sum = 0;
  for (size_t env_id = 0; env_id < ms_map_env_list_.size(); env_id++) {
    MsMapEnvPtr env_ptr = ms_map_env_list_[env_id];
    if (env_ptr == nullptr) {
      continue;
    }
    MsMapPeakPtrVec exp_peak_list = env_ptr->getMsMapPeakList();
    double ratio = ratio_list[env_id];
    for (int peak_idx = 0; peak_idx < num_peaks_theo_env; peak_idx++) {
      MsMapPeakPtr peak = exp_peak_list[peak_idx];
      if (peak != nullptr) {
        double cur_inte = seed_peak_list[peak_idx]->getIntensity() * ratio;
        double cur_err = peak->getPosition() - seed_peak_list[peak_idx]->getPosition();
        error_sum = error_sum + (cur_inte * cur_err);
        weight_sum = weight_sum + cur_inte;
      }
    }
  }
  std::pair<double, double> results(error_sum, weight_sum);
  return results;
}

double getLeftMax(int pos, std::vector<double> &y) {
  double max_val = std::numeric_limits<double>::lowest();
  for (int i = 0; i < pos; i++) {
    if (y[i] > max_val)
      max_val = y[i];
  }
  return max_val;
}

double getRightMax(int pos, std::vector<double> &y) {
  double max_val = std::numeric_limits<double>::lowest();
  int vec_length = y.size();
  for (int i = pos + 1; i < vec_length; i++) {
    if (y[i] > max_val)
      max_val = y[i];
  }
  return max_val;
}

void EnvSet::refineXicBoundary() {
  double split_feature_intensity_ratio = 0.4;
  int seed_idx = seed_ptr_->getSpecId() - start_spec_id_;
  std::vector<double> smoothed_env_xic = xic_ptr_->getSmoothedInteList();

  /// Left side
  std::vector<double> left_data(smoothed_env_xic.begin(), smoothed_env_xic.begin() + seed_idx + 1);
  std::vector<int> minima_left = xic_util::findLocalMinima(left_data);
  std::vector<double> minima_vals_left;
  for (auto m: minima_left) minima_vals_left.push_back(left_data[m]);
  int start_split_point = start_spec_id_;
  while (!minima_vals_left.empty()) {
    int idx = std::min_element(minima_vals_left.begin(), 
                               minima_vals_left.end()) - minima_vals_left.begin();
    int pos = minima_left[idx];
    minima_vals_left.erase(minima_vals_left.begin() + idx);
    double left_max = getLeftMax(pos, left_data);
    if (left_max == 0) continue;
    if (left_data[pos] / left_max <= split_feature_intensity_ratio) {
      start_split_point = start_split_point + pos;
      std::vector<double> temp_left_data(left_data.begin() + pos, left_data.end());
      left_data = temp_left_data;
      minima_left = xic_util::findLocalMinima(left_data);
      minima_vals_left.clear();
      for (auto m: minima_left) minima_vals_left.push_back(left_data[m]);
    }
  }

  /// Right side
  std::vector<double> right_data(smoothed_env_xic.begin() + seed_idx, smoothed_env_xic.end());
  std::vector<int> minima_right = xic_util::findLocalMinima(right_data);
  std::vector<double> minima_vals_right;
  for (auto m: minima_right) minima_vals_right.push_back(right_data[m]);
  int end_split_point = -1;
  while (!minima_vals_right.empty()) {
    int idx = std::min_element(minima_vals_right.begin(), 
                               minima_vals_right.end()) - minima_vals_right.begin();
    int pos = minima_right[idx];
    minima_vals_right.erase(minima_vals_right.begin() + idx);
    double right_max = getRightMax(pos, right_data);
    if (right_max == 0) continue;
    if (right_data[pos] / right_max <= split_feature_intensity_ratio) {
      end_split_point = pos;
      std::vector<double> temp_right_data(right_data.begin(), right_data.begin() + pos - 1);
      right_data = temp_right_data;
      minima_right = xic_util::findLocalMinima(right_data);
      minima_vals_right.clear();
      for (auto m: minima_right) minima_vals_right.push_back(right_data[m]);
    }
  }
  // find new start and end spec ids
  int start = start_spec_id_;
  if (start_split_point > -1)
    start = start_split_point;
  int end = end_spec_id_;
  if (end_split_point > -1)
    end = seed_ptr_->getSpecId() + end_split_point;
  start_spec_id_ = start; 
  end_spec_id_ = end; 

  // update env_set list and xic
  shortlistExpEnvs();
  initMedianXic();
}

bool EnvSet::containTwoValidEnvs(int min_match_peak) {
  int cnt = 0;
  int ref_idx = seed_ptr_->getReferIdx();
  for (size_t i = 0; i < ms_map_env_list_.size(); i++) {
    MsMapEnvPtr env_ptr = ms_map_env_list_[i];
    if (env_ptr != nullptr 
        && env_ptr->getTopThreeMatchNum(ref_idx) >= min_match_peak) {
      cnt++;
    }
  }
  if (cnt >= 2) {
    return true;
  }
  else {
    return false;
  }
}

// check if the seed envelope and one of the neighboring ones are valid
bool EnvSet::containTwoValidOutOfThreeEnvs(int min_match_peak_num) {
  int seed_spec_idx = seed_ptr_->getSpecId() - start_spec_id_; 
  size_t first_idx = seed_spec_idx - 1;
  if (first_idx < 0) {
    first_idx = 0;
  }
  size_t last_idx = seed_spec_idx + 1;
  if (last_idx >= ms_map_env_list_.size()) {
    last_idx = ms_map_env_list_.size() -1;
  }
  int ref_idx = seed_ptr_->getReferIdx(); 
  int count = 0;
  for (size_t i = first_idx; i <= last_idx; i++) {
    MsMapEnvPtr env_ptr = ms_map_env_list_[i];
    if (env_ptr != nullptr && env_ptr->getTopThreeMatchNum(ref_idx) >= min_match_peak_num) {
      count = count + 1;
    }
  }
  if (count >= 2) {
    return true;
  }
  else {
    return false;
  }
}

bool EnvSet::containThreeValidOutOfFiveEnvs(int min_match_peak_num) {
  int seed_spec_idx = seed_ptr_->getSpecId() - start_spec_id_; 
  size_t first_idx = seed_spec_idx - 2;
  if (first_idx < 0) {
    first_idx = 0;
  }
  size_t last_idx = seed_spec_idx + 2;
  if (last_idx >= ms_map_env_list_.size()) {
    last_idx = ms_map_env_list_.size() -1;
  }
  int ref_idx = seed_ptr_->getReferIdx();
  int cnt = 0;
  for (size_t i = first_idx; i <= last_idx; i++) {
    MsMapEnvPtr env_ptr = ms_map_env_list_[i];
    if (env_ptr != nullptr && env_ptr->getTopThreeMatchNum(ref_idx) >= min_match_peak_num) {
      cnt++;
    }
  }
  if (cnt >= 3) {
    return true;
  }
  else {
    return false;
  }
}

void EnvSet::mergeEnvSet(EnvSetPtr new_set_ptr) {
  int new_start_id = new_set_ptr->getStartSpecId();
  int merge_start_id = std::min(start_spec_id_, new_start_id);
  int new_end_id = new_set_ptr->getEndSpecId();
  int merge_end_id = std::max(end_spec_id_, new_end_id);
  // merge 
  MsMapEnvPtrVec merge_env_list (merge_end_id - merge_start_id + 1, nullptr);
  for (size_t i = 0; i < ms_map_env_list_.size(); i++) {
    int idx = start_spec_id_ + i - merge_start_id;
    merge_env_list[idx] = ms_map_env_list_[i];
  }
  MsMapEnvPtrVec new_env_list = new_set_ptr->getMsMapEnvList();
  for (size_t i = 0; i < new_env_list.size(); i++) {
    int idx = new_start_id + i - merge_start_id;
    if (merge_env_list[idx] == nullptr) {
      merge_env_list[idx] = new_env_list[i];
    }
  }
  // assignment
  start_spec_id_ = merge_start_id;
  end_spec_id_ = merge_end_id;
  ms_map_env_list_ = merge_env_list;
  initMedianXic();
}

void EnvSet::appendToXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent) {
  std::string element_name = "env_set";
  XmlDOMElement* element = xml_doc->createElement(element_name.c_str());
  seed_ptr_->appendToXml(xml_doc, element);
  element_name = "exp_env_list";
  XmlDOMElement* env_list = xml_doc->createElement(element_name.c_str());
  for (size_t i = 0; i < ms_map_env_list_.size(); i++) {
    ms_map_env_list_[i]->appendToXml(xml_doc, env_list);
  }
  element->appendChild(env_list);
  parent->appendChild(element);
}

}

