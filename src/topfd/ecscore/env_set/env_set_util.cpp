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
#include <algorithm>

#include "common/util/logger.hpp"

#include "topfd/ecscore/env/seed_env_util.hpp"
#include "topfd/ecscore/env/ms_map_env_util.hpp"
#include "topfd/ecscore/env_set/env_set.hpp"
#include "topfd/ecscore/env_set/env_set_util.hpp"

namespace toppic {

namespace env_set_util {

MsMapPeakPtr pickMsMapPeak(MsMapPtr ms_map_ptr, EnvPeakPtr seed_peak_ptr,
                           int sp_id, double mass_tol) {
  // get peaks within mass tolerance
  double max_inte = std::numeric_limits<double>::min();
  double mz = seed_peak_ptr->getPosition();
  int start_idx = ms_map_ptr->getColIndex(mz - mass_tol);
  if (start_idx < 0) {
      start_idx = 0;
  }
  int end_idx = ms_map_ptr->getColIndex(mz + mass_tol);
  if (end_idx >= ms_map_ptr->getColNum()) {
      end_idx = ms_map_ptr->getColNum() - 1;
  }
  MsMapPeakPtr result_peak = nullptr;
  for (int idx = start_idx; idx <= end_idx; idx++) {
    MsMapPeakPtrVec bin_peaks = ms_map_ptr->getBinPeakList(sp_id, idx);
    for (const auto& peak_ptr : bin_peaks) {
      double mass_diff = std::abs(mz - peak_ptr->getPosition());
      if ( mass_diff < mass_tol && peak_ptr->getIntensity() > max_inte) {
        result_peak = peak_ptr;
        max_inte = peak_ptr->getIntensity();
      }
    }
  }
  return result_peak;
}

MsMapEnvPtr getMatchMsMapEnv(MsMapPtr ms_map_ptr, SeedEnvPtr seed_ptr,
                             int sp_id, double mass_tol) {
  MsMapPeakPtrVec peak_list;
  EnvPeakPtrVec peaks = seed_ptr->getPeakPtrList();
  for (auto& seed_peak : peaks) {
    if (seed_peak != nullptr) {
      MsMapPeakPtr peak = pickMsMapPeak(ms_map_ptr, seed_peak, sp_id, mass_tol);
      peak_list.push_back(peak);
    }
    else {
      peak_list.push_back(nullptr);
    }
  }
  MsMapEnvPtr ms_map_env_ptr = std::make_shared<MsMapEnv>(sp_id, peak_list);
  return ms_map_env_ptr;
}

MsMapEnvPtr getMatchMsMapEnv(MsMapPtr ms_map_ptr, SeedEnvPtr seed_ptr,
                             int sp_id, double mass_tole,
                             std::vector<double> &seed_inte_list,
                             double min_inte) {
  MsMapEnvPtr ms_map_env_ptr = getMatchMsMapEnv(ms_map_ptr, seed_ptr,
                                                sp_id, mass_tole);
  double inte_ratio = ms_map_env_util::compTopThreeInteRatio(seed_ptr, ms_map_env_ptr);
  ms_map_env_ptr->removeLowIntePeaks(seed_inte_list, inte_ratio, min_inte);
  return ms_map_env_ptr;
}

void removeNonMatchEnvs(MsMapEnvPtrVec &env_list, int refer_idx,
                        int min_match_peak_num) {
  int idx = env_list.size() - 1;
  while (idx >= 0) {
    MsMapEnvPtr env = env_list[idx];
    if (env->getTopThreeMatchNum(refer_idx) < min_match_peak_num)
      env_list.erase(env_list.begin() + idx);
    else
      return;
    idx = idx - 1;
  }
}

EnvSetPtr searchEnvSet(MsMapPtr ms_map_ptr, SeedEnvPtr seed_ptr,
                      EcscoreParaPtr para_ptr, double sn_ratio) { 
  int start_spec_id = 0;
  int end_spec_id = ms_map_ptr->getRowNum() - 1;
  return searchEnvSet(ms_map_ptr, seed_ptr, 
                      start_spec_id, end_spec_id,
                      para_ptr, sn_ratio);
}

EnvSetPtr searchEnvSet(MsMapPtr ms_map_ptr, SeedEnvPtr seed_ptr,
                       int start_spec_id, int end_spec_id,
                       EcscoreParaPtr para_ptr, double sn_ratio) {
  double mass_tole = para_ptr->getMassTole();
  int refer_peak_idx = seed_ptr->getReferIdx();
  double min_inte = ms_map_ptr->getBaseInte() * sn_ratio;
  std::vector<double> seed_inte_list = seed_ptr->getInteList();
   // search backward
  MsMapEnvPtrVec back_env_list;
  int miss_num = 0;
  int spec_id = seed_ptr->getSpecId();
  while (spec_id >= start_spec_id) {
    MsMapEnvPtr ms_map_env_ptr = getMatchMsMapEnv(ms_map_ptr, seed_ptr,
                                                  spec_id, mass_tole,
                                                  seed_inte_list, min_inte);
    back_env_list.push_back(ms_map_env_ptr);
    if (ms_map_env_ptr->getTopThreeMatchNum(refer_peak_idx)
       < para_ptr->min_match_peak_) {
      miss_num = miss_num + 1;
    }
    else {
      miss_num = 0;
    }
    if (miss_num >= para_ptr->max_miss_env_) {
      break;
    }
    spec_id = spec_id - 1;
  }
  removeNonMatchEnvs(back_env_list, refer_peak_idx, para_ptr->min_match_peak_);

  // search forward
  MsMapEnvPtrVec forw_env_list;
  spec_id = seed_ptr->getSpecId() + 1;
  miss_num = 0;
  while (spec_id <= end_spec_id) {
    MsMapEnvPtr ms_map_env_ptr = getMatchMsMapEnv(ms_map_ptr, seed_ptr,
                                                  spec_id, mass_tole,
                                                  seed_inte_list, min_inte);
    forw_env_list.push_back(ms_map_env_ptr);
    if (ms_map_env_ptr->getTopThreeMatchNum(refer_peak_idx) < para_ptr->min_match_peak_) {
      miss_num = miss_num + 1;
    }
    else {
      miss_num = 0;
    }
    if (miss_num >= para_ptr->max_miss_env_) {
      break;
    }
    spec_id = spec_id + 1;
  }
  removeNonMatchEnvs(forw_env_list, refer_peak_idx, para_ptr->min_match_peak_);
  // merge results
  std::reverse(back_env_list.begin(), back_env_list.end());
  back_env_list.insert(back_env_list.end(), forw_env_list.begin(), forw_env_list.end());
  if (back_env_list.empty()) {
    return nullptr;
  }
  start_spec_id = back_env_list[0]->getSpecId();
  end_spec_id = back_env_list[back_env_list.size() - 1]->getSpecId();
  if ((end_spec_id - start_spec_id + 1) < para_ptr->min_match_env_) {
    return nullptr;
  }
  EnvSetPtr env_set_ptr = std::make_shared<EnvSet>(seed_ptr, back_env_list, 
                                                   start_spec_id, end_spec_id, 
                                                   min_inte);
  return env_set_ptr;
}

bool checkValidEnvSetSeedEnv(MsMapPtr ms_map_ptr, EnvSetPtr env_set_ptr,
                             int min_match_peak) {
  SeedEnvPtr seed_ptr = env_set_ptr->getSeedPtr();
  std::vector<double> theo_env_inte = seed_ptr->getInteList();
  int refer_idx = std::max_element(theo_env_inte.begin(), theo_env_inte.end()) - theo_env_inte.begin();
  int base_idx = env_set_ptr->getSeedSpecId();
  int start_idx = std::max(base_idx-1, 0);
  int end_idx = std::min(base_idx+1, ms_map_ptr->getRowNum() - 1);
  MsMapEnvPtrVec env_list = env_set_ptr->getMsMapEnvList();
  bool valid = true;
  for (auto &exp_env : env_list) {
    if (exp_env->getSpecId() >= start_idx and exp_env->getSpecId() <= end_idx)
      if (exp_env->getTopThreeMatchNum(refer_idx) < min_match_peak)
        valid = false;
  }
  return valid;
}

bool checkValidEnvSetSeedEnvSparse(MsMapPtr ms_map_ptr, EnvSetPtr env_set_ptr,
                                   int min_match_peak) {
  SeedEnvPtr seed_ptr = env_set_ptr->getSeedPtr();
  std::vector<double> theo_env_inte = seed_ptr->getInteList();
  int refer_idx = std::max_element(theo_env_inte.begin(), theo_env_inte.end()) - theo_env_inte.begin();
  int base_idx = env_set_ptr->getSeedSpecId();
  int start_idx = std::max(base_idx-2, 0);
  int end_idx = std::min(base_idx+2, ms_map_ptr->getRowNum() - 1);
  MsMapEnvPtrVec env_list = env_set_ptr->getMsMapEnvList();
  bool valid = true;
  int false_counter = 0;
  for (auto &exp_env : env_list) {
    if (exp_env->getSpecId() >= start_idx and exp_env->getSpecId() <= end_idx)
      if (exp_env->getTopThreeMatchNum(refer_idx) < min_match_peak)
        false_counter++;
  }
  if (false_counter > 2)
    valid = false;
  return valid;
}

bool checkValidEnvSet(MsMapPtr ms_map_ptr, EnvSetPtr env_set_ptr) {
  bool valid = true;
  int elems = 0;
  std::vector<double> env_xic = env_set_ptr->getXicPtr()->getTopThreeInteList();
  for (double inte : env_xic)
    if (inte > 0) elems++;
  if (elems < 2) valid = false;
  return valid;
}

SeedEnvPtr getHalfChargeEnv(SeedEnvPtr seed_ptr,
                            double even_odd_peak_ratio) {
  double old_charge = seed_ptr->getCharge();
  if (old_charge < 2) {
    return nullptr;
  }
  int new_charge = int(old_charge / 2);
  double ref_mz = seed_ptr->getReferMz();
  double ref_mass = peak_util::compPeakNeutralMass(ref_mz, new_charge);
  double mono_mass = EnvBase::convertRefMassToMonoMass(ref_mass);
  int sp_id = seed_ptr->getSpecId();
  int env_id = -1;
  double inte = seed_ptr->getSeedInte()/2;
  DeconvPeakPtr peak_ptr = std::make_shared<DeconvPeak>(sp_id, env_id,
                                                        mono_mass, inte,
                                                        new_charge);
  SeedEnvPtr new_seed_ptr = std::make_shared<SeedEnv>(peak_ptr);
  return new_seed_ptr;
}

SeedEnvPtr testHalfChargeState(MsMapPtr ms_map_ptr, SeedEnvPtr seed_ptr,
                               EnvSetPtr env_set_ptr, double even_odd_peak_ratio,
                               EcscoreParaPtr para_ptr, double sn_ratio) {
  SeedEnvPtr half_charge_seed = getHalfChargeEnv(seed_ptr, even_odd_peak_ratio);
  bool valid = false;
  valid = seed_env_util::preprocessEnv(ms_map_ptr, half_charge_seed, para_ptr, sn_ratio);
  if (!valid)
    return nullptr;
  return half_charge_seed;
}

}

}
