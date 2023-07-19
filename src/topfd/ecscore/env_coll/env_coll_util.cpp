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
#include <limits>

#include "common/util/logger.hpp"
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
      if (std::count(parent_charge_states.begin(), parent_charge_states.end(),
                     charge_state)) {
        overlap_charge = true;
      }
    }
    if (status) {
      return true;
    }
    /** This part is problematic and needs to be rewritten
    if (status and !overlap_charge) {
      EnvSetPtrVec esl = f->getEnvSetList();
      for (const auto &feature: env_coll_ptr->getEnvSetList()) {
        esl.push_back(feature);
        return true;
      }
      f->setEnvSetList(esl);
    }
    if (status and overlap_charge) {
      return true;
    }
    **/
  }
  return false;
}

EnvSetPtrVec getChargeEnvList(PeakMatrixPtr matrix_ptr, SeedEnvelopePtr seed_ptr,
                              EnvSetPtr seed_env_set_ptr, EcscoreParaPtr para_ptr, 
                              double sn_ratio) {
  int start_spec_id = seed_env_set_ptr->getStartSpecId();
  int end_spec_id = seed_env_set_ptr->getEndSpecId();
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
    if (env_set_ptr == nullptr) {
      miss_num = miss_num + 1;
    } 
    else {
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
    charge = charge - 1;
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

EnvCollPtr findEnvCollWithSingleEnv(PeakMatrixPtr matrix_ptr, SeedEnvelopePtr seed_ptr,
                                    EcscoreParaPtr para_ptr, double sn_ratio) {
  EnvSetPtr env_set_ptr = env_set_util::getEnvSet(matrix_ptr, seed_ptr, para_ptr, sn_ratio);
  if (env_set_ptr == nullptr) {
    LOG_ERROR("env set is null");
    return nullptr; 
  }
  env_set_ptr->refineFeatureBoundary();

  EnvSetPtrVec env_set_list = getChargeEnvList(matrix_ptr, seed_ptr,
                                               env_set_ptr, para_ptr, sn_ratio);
  env_set_list.push_back(env_set_ptr);
  std::sort(env_set_list.begin(), env_set_list.end(), EnvSet::cmpCharge);

  int min_charge = env_set_list[0]->getCharge();
  int max_charge = env_set_list[env_set_list.size() - 1]->getCharge();
  int start_spec_id = env_set_ptr->getStartSpecId();
  int end_spec_id = env_set_ptr->getEndSpecId();
  EnvCollPtr env_coll_ptr = std::make_shared<EnvColl>(seed_ptr, env_set_list, 
                                                      min_charge, max_charge, 
                                                      start_spec_id, end_spec_id);
  return env_coll_ptr;
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
  SeedEnvelopePtr new_seed_ptr = seed_ptr;
  if (std::abs(even_odd_peak_ratio) > para_ptr->even_odd_ratio_cutoff_) {
    new_seed_ptr = env_util::testHalfChargeState(matrix_ptr, seed_ptr,
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
  
  EnvSetPtrVec env_set_list = getChargeEnvList(matrix_ptr, new_seed_ptr,
                                               env_set_ptr, para_ptr, sn_ratio);
  env_set_list.push_back(env_set_ptr);
  std::sort(env_set_list.begin(), env_set_list.end(), EnvSet::cmpCharge);
  int min_charge = env_set_list[0]->getCharge();
  int max_charge = env_set_list[env_set_list.size() - 1]->getCharge();
  int start_spec_id = env_set_ptr->getStartSpecId();
  int end_spec_id = env_set_ptr->getEndSpecId();
  if (new_seed_ptr->getMass() > 10059.3 && new_seed_ptr->getMass() < 10059.4) {
    LOG_ERROR("Mass " << new_seed_ptr->getMass() << " env set list length " <<
              env_set_list.size());
  }

  EnvCollPtr env_coll_ptr = std::make_shared<EnvColl>(new_seed_ptr, env_set_list, 
                                                      min_charge, max_charge, 
                                                      start_spec_id, end_spec_id);
  return env_coll_ptr;
}

FracFeaturePtr getFracFeature(int feat_id, DeconvMsPtrVec &ms1_ptr_vec, int frac_id, 
                              std::string &file_name, EnvCollPtr coll_ptr, 
                              PeakMatrixPtr matrix_ptr, double sn_ratio) {

  double noise_inte = matrix_ptr->getBaseInte(); 
  MatrixSpectrumPtrVec spec_list = matrix_ptr->getSpecList(); 
  int ms1_id_begin = coll_ptr->getStartSpecId();
  int ms1_id_end = coll_ptr->getEndSpecId();
  double feat_inte = coll_ptr->getIntensity(sn_ratio, noise_inte); 
  double feat_mass = coll_ptr->getMass();
  int min_charge = coll_ptr->getMinCharge();
  int max_charge = coll_ptr->getMaxCharge();
  double ms1_time_begin = spec_list[ms1_id_begin]->getRt(); 
  double ms1_time_end = spec_list[ms1_id_end]->getRt(); 
  int ms1_scan_begin = ms1_ptr_vec[ms1_id_begin]->getMsHeaderPtr()->getFirstScanNum();
  int ms1_scan_end = ms1_ptr_vec[ms1_id_end]->getMsHeaderPtr()->getFirstScanNum();
  // get apex inte
  int ms1_apex_id = coll_ptr->getBaseSpecId();
  double time_apex = spec_list[ms1_apex_id]->getRt(); 

  double apex_inte = coll_ptr->getSeedEnvSet()->getXicSeedInte(); 
  FracFeaturePtr feature_ptr = std::make_shared<FracFeature>(feat_id, frac_id, file_name, feat_mass, feat_inte,
                                                             ms1_id_begin, ms1_id_end, ms1_time_begin, ms1_time_end,
                                                             ms1_scan_begin, ms1_scan_end, min_charge, max_charge,
                                                             0, time_apex, apex_inte);
  SingleChargeFeaturePtrVec single_features;
  for (EnvSetPtr es: coll_ptr->getEnvSetList()) {
    int id_begin = es->getStartSpecId();
    int id_end = es->getEndSpecId();
    double time_begin = ms1_ptr_vec[id_begin]->getMsHeaderPtr()->getRetentionTime();
    double time_end = ms1_ptr_vec[id_end]->getMsHeaderPtr()->getRetentionTime();
    int scan_begin = ms1_ptr_vec[id_begin]->getMsHeaderPtr()->getFirstScanNum();
    int scan_end = ms1_ptr_vec[id_end]->getMsHeaderPtr()->getFirstScanNum();
    double inte = es->compIntensity(sn_ratio, noise_inte);
    int env_num = es->countEnvNum();
    int charge = es->getCharge();
    //std::vector<double> xic = es->getXicTopThreeInteList();
    SingleChargeFeaturePtr single_feature = std::make_shared<SingleChargeFeature>(charge, time_begin, time_end,
                                                                                  scan_begin, scan_end,
                                                                                  inte, env_num);
    single_features.push_back(single_feature);
  }
  feature_ptr->setSingleFeatures(single_features);
  return feature_ptr;
}

}
}
