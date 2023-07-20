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

#include <limits>

#include "common/util/logger.hpp"

#include "topfd/ecscore/env_set/env_set.hpp"
#include "topfd/ecscore/env_set/env_set_util.hpp"

namespace toppic {

namespace env_set_util {

std::vector<double> getAggregateEnvelopeInte(EnvSetPtr env_set_ptr) {
  ExpEnvelopePtrVec exp_env_list = env_set_ptr->getExpEnvList();
  std::vector<double> aggregate_inte(exp_env_list[0]->getPeakNum(), 0.0);
  int num_spec = exp_env_list.size();
  int num_peaks = aggregate_inte.size();
  for (int sp_idx = 0; sp_idx < num_spec; sp_idx++) {
    for (int peak_idx = 0; peak_idx < num_peaks; peak_idx++) {
      MsMapPeakPtr peak_ptr = exp_env_list[sp_idx]->getPeakPtr(peak_idx);
      if (peak_ptr != nullptr) {
        aggregate_inte[peak_idx] = aggregate_inte[peak_idx] + peak_ptr->getIntensity();
      }
    }
  }
  return aggregate_inte;
}

std::vector<double> getAggregateEnvelopeMz(EnvSetPtr env_set_ptr) {
  ExpEnvelopePtrVec exp_env_list = env_set_ptr->getExpEnvList();
  std::vector<double> aggregate_mz(exp_env_list[0]->getPeakNum(), 0.0);
  int num_spec = exp_env_list.size();
  int num_peaks = aggregate_mz.size();
  for (int peak_idx = 0; peak_idx < num_peaks; peak_idx++) {
    int counter = 0;
    for (int sp_idx = 0; sp_idx < num_spec; sp_idx++) {
      MsMapPeakPtr peak_ptr = exp_env_list[sp_idx]->getPeakPtr(peak_idx);
      if (peak_ptr != nullptr) {
        aggregate_mz[peak_idx] = aggregate_mz[peak_idx] + peak_ptr->getPosition();
        counter = counter + 1;
      }
    }
    if (counter > 0) {
      aggregate_mz[peak_idx] = aggregate_mz[peak_idx] / counter;
    }
  }
  return aggregate_mz;
}

double calcInteRatio(std::vector<double> &theo_envelope_inte, 
                     std::vector<double> &exp_envelope_inte) {
  if (theo_envelope_inte.size() == 0) {
    LOG_WARN("Empty peak list!");
  }
  double theo_sum = 0;
  double obs_sum = 0;
  int refer_idx =
    std::max_element(theo_envelope_inte.begin(), theo_envelope_inte.end()) - theo_envelope_inte.begin();
  theo_sum = theo_sum + theo_envelope_inte[refer_idx];
  obs_sum = obs_sum + exp_envelope_inte[refer_idx];
  if (refer_idx - 1 >= 0) {
    theo_sum = theo_sum + theo_envelope_inte[refer_idx - 1];
    obs_sum = obs_sum + exp_envelope_inte[refer_idx - 1];
  }
  if (refer_idx + 1 < static_cast<int>(theo_envelope_inte.size())) {
    theo_sum = theo_sum + theo_envelope_inte[refer_idx + 1];
    obs_sum = obs_sum + exp_envelope_inte[refer_idx + 1];
  }
  if (theo_sum == 0) {
    return 1.0;
  }
  else {
    return obs_sum / theo_sum;
  }
}

/*
void compPeakStartEndIdx(MsMapPtr matrix_ptr, SeedEnvelopePtr seed_ptr,
                         double error_tole) {
  SeedEnvPeakPtrVec peak_list = seed_ptr->getPeakList();
  for (auto &peak : peak_list) {
    double mz = peak->getPosition();
    int start_idx = matrix_ptr->getColIndex(mz - error_tole);
    if (start_idx < 0)
      start_idx = 0;
    int end_idx = matrix_ptr->getColIndex(mz + error_tole);
    if (end_idx >= matrix_ptr->getColNum())
      end_idx = matrix_ptr->getColNum() - 1;
    peak->setStartIdx(start_idx);
    peak->setEndIdx(end_idx);
  }
}
*/

MsMapPeakPtr pickExpPeak(MsMapPtr matrix_ptr, EnvPeakPtr seed_peak_ptr,
                         int sp_id, double mass_tol) {
  // get peaks within mass tolerance
  double max_inte = std::numeric_limits<double>::min();
  double mz = seed_peak_ptr->getPosition();
  int start_idx = matrix_ptr->getColIndex(mz - mass_tol);
  if (start_idx < 0) {
      start_idx = 0;
  }
  int end_idx = matrix_ptr->getColIndex(mz + mass_tol);
  if (end_idx >= matrix_ptr->getColNum()) {
      end_idx = matrix_ptr->getColNum() - 1;
  }
  MsMapPeakPtr result_peak = nullptr;
  for (int idx = start_idx; idx <= end_idx; idx++) {
    MsMapPeakPtrVec bin_peaks = matrix_ptr->getBinPeakList(sp_id, idx);
    for (const auto& matrix_peak : bin_peaks) {
      double mass_diff = std::abs(mz - matrix_peak->getPosition());
      if ( mass_diff < mass_tol && matrix_peak->getIntensity() > max_inte) {
        result_peak = matrix_peak; 
        max_inte = matrix_peak->getIntensity();
      }
    }
  }
  return result_peak;
}

ExpEnvelopePtr getMatchExpEnv(MsMapPtr matrix_ptr, SeedEnvelopePtr seed_ptr,
                              int sp_id, double mass_tol) {
  MsMapPeakPtrVec peak_list;
  EnvPeakPtrVec peaks = seed_ptr->getPeakPtrList();
  for (auto& seed_peak : peaks) {
    if (seed_peak != nullptr) {
      MsMapPeakPtr peak = pickExpPeak(matrix_ptr, seed_peak, sp_id, mass_tol);
      peak_list.push_back(peak);
    }
    else {
      peak_list.push_back(nullptr);
    }
  }
  ExpEnvelopePtr exp_env_ptr = std::make_shared<ExpEnvelope>(sp_id, peak_list);
  return exp_env_ptr;
}

void removeNonMatchEnvs(ExpEnvelopePtrVec &env_list, int refer_idx, 
                        int min_match_peak_num) {
  int idx = env_list.size() - 1;
  while (idx >= 0) {
    ExpEnvelopePtr env = env_list[idx];
    if (env->getMatchPeakNum(refer_idx) < min_match_peak_num)
      env_list.erase(env_list.begin() + idx);
    else
      return;
    idx = idx - 1;
  }
}

EnvSetPtr getEnvSet(MsMapPtr matrix_ptr, SeedEnvelopePtr seed_ptr,
                    EcscoreParaPtr para_ptr, double sn_ratio) {
  double noise_inte_level = matrix_ptr->getBaseInte();
  std::vector<double> theo_env_inte = seed_ptr->getInteList();
  int num_theo_env_peaks = theo_env_inte.size();
  int refer_idx = std::max_element(theo_env_inte.begin(), theo_env_inte.end()) - theo_env_inte.begin();
  int idx = seed_ptr->getSpecId();
  // search backward
  ExpEnvelopePtrVec back_env_list;
  int miss_num = 0;
  while (idx >= 0) {
    ExpEnvelopePtr exp_env_ptr = env_set_util::getMatchExpEnv(matrix_ptr, seed_ptr, 
                                                              idx, para_ptr->mass_tole_);
    std::vector<double> exp_env_inte = exp_env_ptr->getInteList();
    double inte_ratio = calcInteRatio(theo_env_inte, exp_env_inte);
    for (int i = 0; i < num_theo_env_peaks; i++) {
      double peak_inte = theo_env_inte[i];
      if ((inte_ratio * peak_inte) < (noise_inte_level * sn_ratio))
        exp_env_ptr->setPeakPtr(i, nullptr);
    }
    back_env_list.push_back(exp_env_ptr);
    if (exp_env_ptr->getMatchPeakNum(refer_idx) < para_ptr->min_match_peak_) {
      miss_num = miss_num + 1;
    }
    else {
      miss_num = 0;
    }
    if (miss_num >= para_ptr->max_miss_env_) {
      break;
    }
    idx = idx - 1;
  }
  removeNonMatchEnvs(back_env_list, refer_idx, para_ptr->min_match_peak_);
  // search forward
  ExpEnvelopePtrVec forw_env_list;
  idx = seed_ptr->getSpecId() + 1;
  miss_num = 0;
  while (idx < matrix_ptr->getRowNum()) {
    ExpEnvelopePtr exp_env_ptr = env_set_util::getMatchExpEnv(matrix_ptr, seed_ptr, 
                                                              idx, para_ptr->mass_tole_);
    std::vector<double> exp_env_inte = exp_env_ptr->getInteList();
    double inte_ratio = calcInteRatio(theo_env_inte, exp_env_inte);
    for (int i = 0; i < num_theo_env_peaks; i++) {
      double peak_inte = theo_env_inte[i];
      if ((inte_ratio * peak_inte) < (noise_inte_level * sn_ratio))
        exp_env_ptr->setPeakPtr(i, nullptr);
    }
    forw_env_list.push_back(exp_env_ptr);
    if (exp_env_ptr->getMatchPeakNum(refer_idx) < para_ptr->min_match_peak_) {
      miss_num = miss_num + 1;
    }
    else {
      miss_num = 0;
    }
    if (miss_num >= para_ptr->max_miss_env_) {
      break;
    }
    idx = idx + 1;
  }
  removeNonMatchEnvs(forw_env_list, refer_idx, para_ptr->min_match_peak_);
  // merge results
  std::reverse(back_env_list.begin(), back_env_list.end());
  back_env_list.insert(back_env_list.end(), forw_env_list.begin(), forw_env_list.end());
  if (back_env_list.empty()) {
    return nullptr;
  }
  int start_spec_id = back_env_list[0]->getSpecId();
  int end_spec_id = back_env_list[back_env_list.size() - 1]->getSpecId();
  if ((end_spec_id - start_spec_id + 1) < para_ptr->min_match_env_) {
    return nullptr;
  }
  EnvSetPtr env_set_ptr = std::make_shared<EnvSet>(seed_ptr, back_env_list, 
                                                   start_spec_id, end_spec_id, 
                                                   noise_inte_level, sn_ratio);
  return env_set_ptr;
}

bool checkValidEnvSetSeedEnv(MsMapPtr matrix_ptr, EnvSetPtr env_set_ptr,
                             int min_match_peak) {
  std::vector<double> theo_env_inte = env_set_ptr->getSeedInteList();
  int refer_idx = std::max_element(theo_env_inte.begin(), theo_env_inte.end()) - theo_env_inte.begin();
  int base_idx = env_set_ptr->getBaseSpecId();
  int start_idx = std::max(base_idx-1, 0);
  int end_idx = std::min(base_idx+1, matrix_ptr->getRowNum() - 1);
  ExpEnvelopePtrVec env_list = env_set_ptr->getExpEnvList();
  bool valid = true;
  for (auto &exp_env : env_list) {
    if (exp_env->getSpecId() >= start_idx and exp_env->getSpecId() <= end_idx)
      if (exp_env->getMatchPeakNum(refer_idx) < min_match_peak)
        valid = false;
  }
  return valid;
}

bool checkValidEnvSetSeedEnvSparse(MsMapPtr matrix_ptr, EnvSetPtr env_set_ptr,
                                   int min_match_peak) {
  std::vector<double> theo_env_inte = env_set_ptr->getSeedInteList();
  int refer_idx = std::max_element(theo_env_inte.begin(), theo_env_inte.end()) - theo_env_inte.begin();
  int base_idx = env_set_ptr->getBaseSpecId();
  int start_idx = std::max(base_idx-2, 0);
  int end_idx = std::min(base_idx+2, matrix_ptr->getRowNum() - 1);
  ExpEnvelopePtrVec env_list = env_set_ptr->getExpEnvList();
  bool valid = true;
  int false_counter = 0;
  for (auto &exp_env : env_list) {
    if (exp_env->getSpecId() >= start_idx and exp_env->getSpecId() <= end_idx)
      if (exp_env->getMatchPeakNum(refer_idx) < min_match_peak)
        false_counter++;
  }
  if (false_counter > 2)
    valid = false;
  return valid;
}




EnvSetPtr findEnvSet(MsMapPtr matrix_ptr, SeedEnvelopePtr seed_ptr,
                     int start_spec_id, int end_spec_id,
                     EcscoreParaPtr para_ptr, double sn_ratio) {

  double noise_inte_level = matrix_ptr->getBaseInte();
  std::vector<double> theo_envelope_inte = seed_ptr->getInteList();
  int num_theo_env_peaks = theo_envelope_inte.size();
  int refer_idx = std::max_element(theo_envelope_inte.begin(), theo_envelope_inte.end()) - theo_envelope_inte.begin();
  int base_idx = seed_ptr->getSpecId();
  int miss_num = 0;

  ExpEnvelopePtrVec back_env_list;
  for (int idx = base_idx; idx >= start_spec_id; idx--) {
    ExpEnvelopePtr exp_env = env_set_util::getMatchExpEnv(matrix_ptr, seed_ptr, 
                                                          idx, para_ptr->mass_tole_);
    std::vector<double> experimental_envelope_inte = exp_env->getInteList();
    double inte_ratio = calcInteRatio(theo_envelope_inte, experimental_envelope_inte);
    for (int i = 0; i < num_theo_env_peaks; i++) {
      double peak_inte = theo_envelope_inte[i];
      if ((inte_ratio * peak_inte) < (noise_inte_level * sn_ratio))
        exp_env->setPeakPtr(i, nullptr); 
    }
    back_env_list.push_back(exp_env);
    if (exp_env->getMatchPeakNum(refer_idx) < para_ptr->max_miss_peak_)
      miss_num = miss_num + 1;
    else
      miss_num = 0;
    if (miss_num >=  para_ptr->max_miss_env_)
      break;
  }
  removeNonMatchEnvs(back_env_list, refer_idx, para_ptr->min_match_peak_);

  ExpEnvelopePtrVec forw_env_list;
  for (int idx = base_idx + 1; idx <= end_spec_id; idx++) {
    ExpEnvelopePtr exp_env = env_set_util::getMatchExpEnv(matrix_ptr, seed_ptr, idx,  para_ptr->mass_tole_);
    std::vector<double> experimental_envelope_inte = exp_env->getInteList();
    double inte_ratio = calcInteRatio(theo_envelope_inte, experimental_envelope_inte);
    for (int i = 0; i < num_theo_env_peaks; i++) {
      double peak_inte = theo_envelope_inte[i];
      if ((inte_ratio * peak_inte) < (noise_inte_level * sn_ratio))
        exp_env->setPeakPtr(i, nullptr); 
    }
    forw_env_list.push_back(exp_env);
    if (exp_env->getMatchPeakNum(refer_idx) <  para_ptr->max_miss_peak_)
      miss_num = miss_num + 1;
    else
      miss_num = 0;
    if (miss_num >= para_ptr->max_miss_env_)
      break;
  }
  removeNonMatchEnvs(forw_env_list, refer_idx, para_ptr->min_match_peak_);
  // merge
  std::reverse(back_env_list.begin(), back_env_list.end());
  back_env_list.insert(back_env_list.end(), forw_env_list.begin(), forw_env_list.end());
  if (back_env_list.empty()) return nullptr;
  start_spec_id = back_env_list[0]->getSpecId();
  end_spec_id = back_env_list[back_env_list.size() - 1]->getSpecId();
  if ((end_spec_id - start_spec_id + 1) < 2) return nullptr;
  EnvSetPtr env_set_ptr = std::make_shared<EnvSet>(seed_ptr, back_env_list, 
                                                   start_spec_id, end_spec_id, 
                                                   noise_inte_level, sn_ratio);
  return env_set_ptr;
}

bool checkValidEnvSet(MsMapPtr matrix_ptr, EnvSetPtr env_set_ptr) {
  bool valid = true;
  int elems = 0;
  std::vector<double> env_xic = env_set_ptr->getXicTopThreeInteList();
  for (double inte : env_xic)
    if (inte > 0) elems++;
  if (elems < 2) valid = false;
  return valid;
}

}

}



