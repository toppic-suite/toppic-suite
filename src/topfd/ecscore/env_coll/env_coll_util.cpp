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

#include "common/base/mass_constant.hpp"
#include "topfd/ecscore/envelope/env_util.hpp"
#include "topfd/ecscore/env_set/env_set_util.hpp"
#include "topfd/ecscore/score/component_score.hpp"
#include "topfd/ecscore/env_coll/env_coll_util.hpp"

namespace toppic {
namespace env_coll_util {

void mergeOverlappingChargeFeatures(EnvCollPtr f, EnvCollPtr env_coll) {
  std::vector<int> parent_charge_states = f->getChargeList();
  EnvSetPtrVec esl = f->getEnvSetList();
  for (auto &feature: env_coll->getEnvSetList())
    if (std::count(parent_charge_states.begin(), parent_charge_states.end(),
                   feature->getCharge()))
      esl.push_back(feature);
  f->setEnvSetList(esl);
}

bool checkChargeStateDist(const std::vector<int> &parent_charge_states, 
                          int charge_state) {
  int min_charge_diff = std::numeric_limits<int>::max();
  for (auto parent_charge_state: parent_charge_states) {
    int charge_diff = std::abs(charge_state - parent_charge_state);
    if (charge_diff < min_charge_diff)
      min_charge_diff = charge_diff;
  }
  if (min_charge_diff > 2) {
    return false;
  }
  else {
    return true;
  }
}


bool checkOverlap(MatrixSpectrumPtrVec &spectrum_list, EnvCollPtr f, 
                  double feature_start_rt, 
                  double feature_end_rt, double time_tol) {
  int start_spec_id = f->getStartSpecId();
  int end_spec_id = f->getEndSpecId();
  double start_rt = spectrum_list[start_spec_id]->getRt();
  double end_rt = spectrum_list[end_spec_id]->getRt();
  double start = -1;
  if (start_rt <= feature_start_rt)
    start = feature_start_rt;
  else
    start = start_rt;
  double end = -1;
  if (start > -1) {
    if (end_rt <= feature_end_rt)
      end = end_rt;
    else
      end = feature_end_rt;
  }
  bool status = false;
  if (end > -1) {
    double overlapping_rt_range = end - start;
    if (overlapping_rt_range > 0) {
      double feature_rt_range = feature_end_rt - feature_start_rt;
      double feature_coverage = overlapping_rt_range / feature_rt_range;
      if (feature_coverage > time_tol)
        status = true;
    }
  }
  return status;
}

bool checkExistingFeatures(PeakMatrixPtr matrix_ptr, EnvCollPtr env_coll_ptr, 
                           EnvCollPtrVec &env_coll_list, EcscoreParaPtr para_ptr) {
  double env_mass = env_coll_ptr->getMass();
  double mass_tol = para_ptr->match_envelope_tolerance_ * env_mass;
  std::vector<int> charge_states = env_coll_ptr->getChargeList();
  double isotope_mass = mass_constant::getIsotopeMass();
  std::vector<double> extended_masses = {env_mass - isotope_mass, env_mass, env_mass + isotope_mass};
  int num_env_colls = env_coll_list.size();
  MatrixSpectrumPtrVec spectrum_list = matrix_ptr->getSpecList();
  int start_spec_id = env_coll_ptr->getStartSpecId();
  int end_spec_id = env_coll_ptr->getEndSpecId();
  double feature_start_rt = spectrum_list[start_spec_id]->getRt();
  double feature_end_rt = spectrum_list[end_spec_id]->getRt();

  std::vector<int> selected_features;
  for (int i = 0; i < num_env_colls; i++) {
    bool overlap = checkOverlap(spectrum_list, env_coll_list[i], feature_start_rt, feature_end_rt, 
                                para_ptr->time_overlap_tole_); 
    if (overlap) {
      double min_mass_diff = std::numeric_limits<double>::max();
      for (auto ext_mass: extended_masses) {
        double mass_diff = std::abs(ext_mass - env_coll_list[i]->getMass());
        if (mass_diff < min_mass_diff)
          min_mass_diff = mass_diff;
      }
      if (min_mass_diff < mass_tol)
        selected_features.push_back(i);
    }
  }
  bool status = true;
  bool overlap_charge = false;
  for (auto f_idx: selected_features) {
    EnvCollPtr f = env_coll_list[f_idx];
    std::vector<int> parent_charge_states = f->getChargeList();
    for (auto charge_state: charge_states) {
      status = checkChargeStateDist(parent_charge_states, charge_state);
      if (std::count(parent_charge_states.begin(), parent_charge_states.end(), charge_state))
        overlap_charge = true;
    }
    if (status and !overlap_charge) {
      EnvSetPtrVec esl = f->getEnvSetList();
      for (const auto &feature: env_coll_ptr->getEnvSetList()) {
        esl.push_back(feature);
        return true;
      }
      f->setEnvSetList(esl);
    }
    if (status and overlap_charge) {
      mergeOverlappingChargeFeatures(f, env_coll_ptr);
      return true;
    }
  }
  return false;
}

EnvSetPtrVec getChargeEnvList(PeakMatrixPtr matrix_ptr, SeedEnvelopePtr seed_ptr,
                              EnvSetPtr env_set_ptr, EcscoreParaPtr para_ptr, 
                              double sn_ratio) {
  int start_spec_id = env_set_ptr->getStartSpecId();
  int end_spec_id = env_set_ptr->getEndSpecId();
  EnvSetPtrVec env_set_list;
  int charge = seed_ptr->getCharge() - 1;
  int miss_num = 0;
  while (charge >= para_ptr->para_min_charge_) {
    SeedEnvelopePtr cur_seed_ptr = std::make_shared<SeedEnvelope>(seed_ptr, charge); 
    env_set_util::compPeakStartEndIdx(matrix_ptr, cur_seed_ptr, 
                                      para_ptr->mass_tole_);
    EnvSetPtr env_set_ptr = env_set_util::findEnvSet(matrix_ptr, cur_seed_ptr, 
                                                     start_spec_id, end_spec_id, 
                                                     para_ptr, sn_ratio);
    charge = charge - 1;
    if (env_set_ptr == nullptr) {
      miss_num = miss_num + 1;
    } else {
      env_set_ptr->refineFeatureBoundary();
      if (!env_set_util::checkValidEnvSet(matrix_ptr, env_set_ptr))
        miss_num = miss_num + 1;
      else {
        miss_num = 0;
        env_set_list.push_back(env_set_ptr);
      }
    }
    if (miss_num >= para_ptr->max_miss_charge_)
      break;
  }

  miss_num = 0;
  charge = seed_ptr->getCharge() + 1;
  while (charge <= para_ptr->para_max_charge_) {
    SeedEnvelopePtr cur_seed_ptr = std::make_shared<SeedEnvelope>(seed_ptr, charge); 
    env_set_util::compPeakStartEndIdx(matrix_ptr, cur_seed_ptr, para_ptr->mass_tole_);
    EnvSetPtr env_set_ptr = env_set_util::findEnvSet(matrix_ptr, cur_seed_ptr, 
                                                     start_spec_id, end_spec_id, 
                                                     para_ptr, sn_ratio);
    charge = charge + 1;
    if (env_set_ptr == nullptr) {
      miss_num = miss_num + 1;
    } else {
      env_set_ptr->refineFeatureBoundary();
      if (!env_set_util::checkValidEnvSet(matrix_ptr, env_set_ptr))
        miss_num = miss_num + 1;
      else {
        miss_num = 0;
        env_set_list.push_back(env_set_ptr);
      }
    }
    if (miss_num >= para_ptr->max_miss_charge_)
      break;
  }
  std::sort(env_set_list.begin(), env_set_list.end(), EnvSet::cmpCharge);
  return env_set_list;
}

EnvCollPtr findEnvColl(PeakMatrixPtr matrix_ptr, SeedEnvelopePtr seed_ptr,
                       EcscoreParaPtr para_ptr, double sn_ratio) {
  EnvSetPtr env_set_ptr = env_set_util::getEnvSet(matrix_ptr, seed_ptr, para_ptr, sn_ratio);

  if (env_set_ptr == nullptr) {
    return nullptr; 
  }
  env_set_ptr->refineFeatureBoundary();
  if (para_ptr->filter_neighboring_peaks_) {
    if (!env_set_util::checkValidEnvSetSeedEnv(matrix_ptr, env_set_ptr, 
                                               para_ptr->max_miss_peak_))
      return nullptr; 
    else
      if (!env_set_util::checkValidEnvSetSeedEnvSparse(matrix_ptr, env_set_ptr,
                                                       para_ptr->max_miss_peak_))
        return nullptr;
  }
  double even_odd_peak_ratio = component_score::getAggOddEvenPeakRatio(env_set_ptr);
  if (std::abs(even_odd_peak_ratio) > para_ptr->even_odd_ratio_cutoff_) {
    SeedEnvelopePtr new_seed_ptr = env_util::testHalfChargeState(matrix_ptr, seed_ptr,
                                                                 env_set_ptr, even_odd_peak_ratio, 
                                                                 para_ptr, sn_ratio);
    if (new_seed_ptr == nullptr) {
      return nullptr;
    }
    EnvSetPtr tmp_env_set_ptr = env_set_util::getEnvSet(matrix_ptr, new_seed_ptr, 
                                                        para_ptr, sn_ratio);
    if (tmp_env_set_ptr == nullptr) {
      return nullptr;
    }
    env_set_ptr = tmp_env_set_ptr; 
    env_set_ptr->refineFeatureBoundary();
    if (!env_set_util::checkValidEnvSetSeedEnv(matrix_ptr, env_set_ptr,
                                               para_ptr->max_miss_peak_)) {
      return nullptr; 
    }
  }
  EnvSetPtrVec env_set_list = getChargeEnvList(matrix_ptr, seed_ptr,
                                               env_set_ptr, para_ptr, sn_ratio);
  env_set_list.push_back(env_set_ptr);
  std::sort(env_set_list.begin(), env_set_list.end(), EnvSet::cmpCharge);
  int min_charge = env_set_list[0]->getCharge();
  int max_charge = env_set_list[env_set_list.size() - 1]->getCharge();
  if (env_set_list.empty()) {
    return nullptr;
  }
  int start_spec_id = env_set_ptr->getStartSpecId();
  int end_spec_id = env_set_ptr->getEndSpecId();
  EnvCollPtr env_coll_ptr = std::make_shared<EnvColl>(seed_ptr, env_set_list, 
                                                      min_charge, max_charge, 
                                                      start_spec_id, end_spec_id);
  return env_coll_ptr;
}

}
}
