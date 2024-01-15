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
#include "common/base/mass_constant.hpp"
#include "topfd/ecscore/env/seed_env_util.hpp"
#include "topfd/ecscore/env/ms_map_env_util.hpp"
#include "topfd/ecscore/env_set/env_set_util.hpp"
#include "topfd/ecscore/score/component_score.hpp"
#include "topfd/ecscore/env_coll/env_coll_util.hpp"

namespace toppic {

namespace env_coll_util {

EnvCollPtr getEnvCollPtr(MsMapPtr matrix_ptr, SeedEnvPtr seed_ptr,
                         EnvSetPtr seed_env_set_ptr, EcscoreParaPtr para_ptr,
                         double sn_ratio) {
  int start_row_id = seed_env_set_ptr->getStartRowId();
  int end_row_id = seed_env_set_ptr->getEndRowId();
  EnvSetPtrVec env_set_list;
  int charge = seed_ptr->getCharge() - 1;
  int miss_num = 0;
  int min_match_peak_num = para_ptr->getMinMatchPeakNumInTopThree();
  int min_scan_num = para_ptr->min_scan_num_; 
  while (charge >= para_ptr->para_min_charge_) {
    SeedEnvPtr cur_seed_ptr = std::make_shared<SeedEnv>(seed_ptr, charge);
    EnvSetPtr env_set_ptr = env_set_util::searchEnvSet(matrix_ptr, cur_seed_ptr,
                                                       start_row_id, end_row_id,
                                                       para_ptr, sn_ratio);
    if (env_set_ptr == nullptr) {
      miss_num = miss_num + 1;
    } 
    else {
      env_set_ptr->refineXicBoundary();
      if (!env_set_ptr->containValidEnvs(min_scan_num, min_match_peak_num)) {
        miss_num = miss_num + 1;
      }
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
    SeedEnvPtr cur_seed_ptr = std::make_shared<SeedEnv>(seed_ptr, charge);
    EnvSetPtr env_set_ptr = env_set_util::searchEnvSet(matrix_ptr, cur_seed_ptr,
                                                       start_row_id, end_row_id,
                                                       para_ptr, sn_ratio);
    charge = charge + 1;
    if (env_set_ptr == nullptr) {
      miss_num = miss_num + 1;
    } else {
      env_set_ptr->refineXicBoundary();
      if (!env_set_ptr->containValidEnvs(min_scan_num, min_match_peak_num)) {
        miss_num = miss_num + 1;
      }
      else {
        miss_num = 0;
        env_set_list.push_back(env_set_ptr);
      }
    }
    if (miss_num >= para_ptr->max_miss_charge_)
      break;
  }
  env_set_list.push_back(seed_env_set_ptr);
  std::sort(env_set_list.begin(), env_set_list.end(), EnvSet::cmpChargeInc);
  int min_charge = env_set_list[0]->getCharge();
  int max_charge = env_set_list[env_set_list.size() - 1]->getCharge();
  EnvCollPtr env_coll_ptr = std::make_shared<EnvColl>(env_set_list, seed_ptr,  
                                                      min_charge, max_charge,
                                                      start_row_id, end_row_id);

  return env_coll_ptr;
}

EnvCollPtr findEnvColl(MsMapPtr matrix_ptr, SeedEnvPtr seed_ptr,
                       EcscoreParaPtr para_ptr, double sn_ratio) {
  EnvSetPtr env_set_ptr = env_set_util::searchEnvSet(matrix_ptr, seed_ptr, para_ptr, sn_ratio);
  if (env_set_ptr == nullptr) {
    return nullptr; 
  }
  env_set_ptr->refineXicBoundary();
  int min_match_peak_num = para_ptr->getMinMatchPeakNumInTopThree();
  int min_scan_num = para_ptr->min_scan_num_; 
  // check if there are valid envelopes 
  if (!env_set_ptr->containValidEnvs(min_scan_num, min_match_peak_num)) {
    return nullptr;
  }
  double even_odd_peak_ratio = component_score::getAggOddEvenPeakRatio(env_set_ptr);
  SeedEnvPtr new_seed_ptr = seed_ptr;
  if (std::abs(even_odd_peak_ratio) > para_ptr->even_odd_ratio_cutoff_) {
    new_seed_ptr = seed_env_util::testHalfChargeEnv(seed_ptr, matrix_ptr, 
                                                    even_odd_peak_ratio, 
                                                    para_ptr, sn_ratio);
    if (new_seed_ptr == nullptr) {
      return nullptr;
    }
    EnvSetPtr tmp_env_set_ptr = env_set_util::searchEnvSet(matrix_ptr, new_seed_ptr,
                                                           para_ptr, sn_ratio);
    if (tmp_env_set_ptr == nullptr) {
      return nullptr;
    }
    env_set_ptr = tmp_env_set_ptr;
    env_set_ptr->refineXicBoundary();

    if (!env_set_ptr->containValidEnvs(min_scan_num, min_match_peak_num)) {
      return nullptr; 
    }
  }
  
  EnvCollPtr env_coll_ptr = getEnvCollPtr(matrix_ptr, new_seed_ptr, 
                                          env_set_ptr, para_ptr, sn_ratio);
  
  return env_coll_ptr;
}

EnvCollPtr findEnvCollWithSingleEnv(MsMapPtr matrix_ptr, SeedEnvPtr seed_ptr,
                                    EcscoreParaPtr para_ptr, double sn_ratio) {
  EnvSetPtr env_set_ptr = env_set_util::searchEnvSet(matrix_ptr, seed_ptr, para_ptr, sn_ratio);
  if (env_set_ptr == nullptr) {
    LOG_INFO("Envelope set is empty!");
    return nullptr; 
  }
  env_set_ptr->refineXicBoundary();

  EnvCollPtr env_coll_ptr = getEnvCollPtr(matrix_ptr, seed_ptr, 
                                          env_set_ptr, para_ptr, sn_ratio);
  return env_coll_ptr;
}

bool checkOverlap(MsMapRowHeaderPtrVec &spectrum_list, EnvCollPtr coll_ptr,
                  double feature_start_rt,
                  double feature_end_rt, double time_tol) {
  int start_row_id = coll_ptr->getStartRowId();
  int end_row_id = coll_ptr->getEndRowId();
  double start_rt = spectrum_list[start_row_id]->getRt();
  double end_rt = spectrum_list[end_row_id]->getRt();
  double start = -1;
  if (start_rt <= feature_start_rt) {
    start = feature_start_rt;
  }
  else {
    start = start_rt;
  }
  double end = -1;
  if (start > -1) {
    if (end_rt <= feature_end_rt) {
      end = end_rt;
    }
    else {
      end = feature_end_rt;
    }
  }
  if (end > -1) {
    // if the coverage is 100%. Sometimes start = end for single scan
    // and the coverage is 100%.
    if (start == feature_start_rt && end == feature_end_rt) {
      return true;
    }
    double overlapping_rt_range = end - start;
    if (overlapping_rt_range > 0) {
      double feature_rt_range = feature_end_rt - feature_start_rt;
      double feature_coverage = overlapping_rt_range / feature_rt_range;
      // return overlap
      if (feature_coverage > time_tol) {
        return true;
      }
    }
  }
  return false;
}

bool checkExistingFeatures(MsMapPtr matrix_ptr, EnvCollPtr env_coll_ptr,
                           EnvCollPtrVec &env_coll_list, EcscoreParaPtr para_ptr) {
  double env_mass = env_coll_ptr->getMonoNeutralMass();
  double mass_tol = para_ptr->match_feature_ppm_tolerance_ * env_mass;
  std::vector<int> charge_states = env_coll_ptr->getChargeList();
  double isotope_mass = mass_constant::getIsotopeMass();
  std::vector<double> extended_masses; 
  if (env_mass > 5000) {
    extended_masses = {env_mass - isotope_mass, env_mass, env_mass + isotope_mass};
  }
  else {
    extended_masses = {env_mass};
  }
  int num_env_colls = env_coll_list.size();
  MsMapRowHeaderPtrVec spectrum_list = matrix_ptr->getHeaderPtrList();
  int start_row_id = env_coll_ptr->getStartRowId();
  int end_row_id = env_coll_ptr->getEndRowId();
  double feature_start_rt = spectrum_list[start_row_id]->getRt();
  double feature_end_rt = spectrum_list[end_row_id]->getRt();

  EnvCollPtr overlap_env_coll_ptr;
  for (int i = 0; i < num_env_colls; i++) {
    bool overlap = checkOverlap(spectrum_list, env_coll_list[i], feature_start_rt, feature_end_rt, 
                                para_ptr->match_feature_time_overlap_tole_); 
    if (overlap) {
      double min_mass_diff = std::numeric_limits<double>::max();
      for (auto ext_mass: extended_masses) {
        double mass_diff = std::abs(ext_mass - env_coll_list[i]->getMonoNeutralMass());
        if (mass_diff < min_mass_diff)
          min_mass_diff = mass_diff;
      }
      if (min_mass_diff < mass_tol) {
        overlap_env_coll_ptr = env_coll_list[i];
        break;
      }
    }
  }

  if (overlap_env_coll_ptr != nullptr) {
    //EnvSetPtrVec merge_set_ptrs = overlap_env_coll_ptr->getEnvSetList();
    EnvSetPtrVec new_set_ptrs = env_coll_ptr->getEnvSetList();
    //merge_set_ptrs.insert(merge_set_ptrs.end(), new_set_ptrs.begin(), new_set_ptrs.end());
    //overlap_env_coll_ptr->setEnvSetList(merge_set_ptrs);
    for (size_t i = 0; i < new_set_ptrs.size(); i++) {
      overlap_env_coll_ptr->mergeEnvSet(new_set_ptrs[i]);
    }
    return true;
  }
  else {
    return false;
  }
}

bool checkSeedExistingFeatures(MsMapPtr matrix_ptr, SeedEnvPtr seed_ptr,
                               EnvCollPtrVec &env_coll_list, EcscoreParaPtr para_ptr) {
  double env_mass = seed_ptr->getMonoNeutralMass();
  double mass_tol = para_ptr->match_feature_ppm_tolerance_ * env_mass;
  double isotope_mass = mass_constant::getIsotopeMass();
  std::vector<double> extended_masses; 
  if (env_mass > 5000) {
    extended_masses = {env_mass - isotope_mass, env_mass, env_mass + isotope_mass};
  }
  else {
    extended_masses = {env_mass};
  }
  int row_id = seed_ptr->getRowId(); 
  MsMapRowHeaderPtrVec spectrum_list = matrix_ptr->getHeaderPtrList();
  double rt = spectrum_list[row_id]->getRt();
  int num_env_colls = env_coll_list.size();
  EnvCollPtr overlap_env_coll_ptr;
  for (int i = 0; i < num_env_colls; i++) {
    bool overlap = checkOverlap(spectrum_list, env_coll_list[i], rt, rt, 
                                para_ptr->match_feature_time_overlap_tole_); 
    if (overlap) {
      double min_mass_diff = std::numeric_limits<double>::max();
      for (auto ext_mass: extended_masses) {
        double mass_diff = std::abs(ext_mass - env_coll_list[i]->getMonoNeutralMass());
        if (mass_diff < min_mass_diff)
          min_mass_diff = mass_diff;
      }
      if (min_mass_diff < mass_tol) {
        overlap_env_coll_ptr = env_coll_list[i];
        break;
      }
    }
  }

  if (overlap_env_coll_ptr != nullptr) {
    return true;
  }
  else {
    return false;
  }
}

FracFeaturePtr getFracFeature(int feat_id, DeconvMsPtrVec &ms1_ptr_vec, int frac_id,
                              std::string &file_name, EnvCollPtr coll_ptr,
                              MsMapPtr matrix_ptr, double sn_ratio) {

  MsMapRowHeaderPtrVec spec_list = matrix_ptr->getHeaderPtrList();
  int ms1_row_begin = coll_ptr->getStartRowId();
  int ms1_row_end = coll_ptr->getEndRowId();
  double feat_inte = coll_ptr->getIntensity();
  double feat_mass = coll_ptr->getMonoNeutralMass();
  int min_charge = coll_ptr->getMinCharge();
  int max_charge = coll_ptr->getMaxCharge();
  int ms1_id_begin = spec_list[ms1_row_begin]->getRawSpecId(); 
  double ms1_id_end = spec_list[ms1_row_end]->getRawSpecId(); 
  double ms1_time_begin = spec_list[ms1_row_begin]->getRt(); 
  double ms1_time_end = spec_list[ms1_row_end]->getRt(); 
  int ms1_scan_begin = ms1_ptr_vec[ms1_row_begin]->getMsHeaderPtr()->getFirstScanNum();
  int ms1_scan_end = ms1_ptr_vec[ms1_row_end]->getMsHeaderPtr()->getFirstScanNum();
  // get apex inte
  int ms1_apex_row = coll_ptr->getSeedRowId();
  double apex_time = spec_list[ms1_apex_row]->getRt(); 
  int apex_scan = spec_list[ms1_apex_row]->getScanNum(); 
  double apex_inte = coll_ptr->getSeedEnvSet()->getXicSeedAllPeakInte();

  int rep_charge = coll_ptr->getSeedPtr()->getCharge(); 
  double rep_avg_mz = coll_ptr->getSeedPtr()->getAvgMz(); 
  int env_num = coll_ptr->countEnvNum();
  double ec_score = coll_ptr->getEcscore();

  FracFeaturePtr feature_ptr = std::make_shared<FracFeature>(feat_id, frac_id, file_name, feat_mass, feat_inte,
                                                             ms1_id_begin, ms1_id_end, ms1_time_begin, ms1_time_end,
                                                             ms1_scan_begin, ms1_scan_end, min_charge, max_charge,
                                                             apex_time, apex_scan, apex_inte, rep_charge, rep_avg_mz, 
                                                             env_num, ec_score);
  SingleChargeFeaturePtrVec single_features;
  for (EnvSetPtr es: coll_ptr->getEnvSetList()) {
    int row_begin = es->getStartRowId();
    int row_end = es->getEndRowId();
    double time_begin = ms1_ptr_vec[row_begin]->getMsHeaderPtr()->getRetentionTime();
    double time_end = ms1_ptr_vec[row_end]->getMsHeaderPtr()->getRetentionTime();
    int scan_begin = ms1_ptr_vec[row_begin]->getMsHeaderPtr()->getFirstScanNum();
    int scan_end = ms1_ptr_vec[row_end]->getMsHeaderPtr()->getFirstScanNum();
    double inte = es->getInte();
    int env_num = es->countEnvNum();
    int charge = es->getCharge();
    SingleChargeFeaturePtr single_feature = std::make_shared<SingleChargeFeature>(charge, time_begin, time_end,
                                                                                  scan_begin, scan_end,
                                                                                  inte, env_num);
    single_features.push_back(single_feature);
  }
  feature_ptr->setSingleFeatures(single_features);
  return feature_ptr;
}

FracFeaturePtr getFracFeature(int feat_id, DeconvMsPtrVec &ms1_ptr_vec, int frac_id,
                              std::string &file_name, SeedEnvPtr seed_ptr, 
                              MsMapPtr matrix_ptr) {

  MsMapRowHeaderPtrVec spec_list = matrix_ptr->getHeaderPtrList();
  int ms1_row_begin = seed_ptr->getRowId();
  int ms1_row_end = seed_ptr->getRowId();
  double feat_inte = seed_ptr->compInteSum();
  double feat_mass = seed_ptr->getMonoNeutralMass();
  int min_charge = seed_ptr->getCharge();
  int max_charge = seed_ptr->getCharge();
  int ms1_id_begin = spec_list[ms1_row_begin]->getRawSpecId(); 
  double ms1_id_end = spec_list[ms1_row_end]->getRawSpecId(); 
  double ms1_time_begin = spec_list[ms1_row_begin]->getRt(); 
  double ms1_time_end = spec_list[ms1_row_end]->getRt(); 
  int ms1_scan_begin = ms1_ptr_vec[ms1_row_begin]->getMsHeaderPtr()->getFirstScanNum();
  int ms1_scan_end = ms1_ptr_vec[ms1_row_end]->getMsHeaderPtr()->getFirstScanNum();
  // get apex inte
  int ms1_apex_row = seed_ptr->getRowId();
  double apex_time = spec_list[ms1_apex_row]->getRt(); 
  int apex_scan = spec_list[ms1_apex_row]->getScanNum(); 
  double apex_inte = feat_inte;

  int rep_charge = seed_ptr->getCharge(); 
  double rep_avg_mz = seed_ptr->getAvgMz(); 
  int env_num = 1; 
  double ec_score = 0;

  FracFeaturePtr feature_ptr = std::make_shared<FracFeature>(feat_id, frac_id, file_name, feat_mass, feat_inte,
                                                             ms1_id_begin, ms1_id_end, ms1_time_begin, ms1_time_end,
                                                             ms1_scan_begin, ms1_scan_end, min_charge, max_charge,
                                                             apex_time, apex_scan, apex_inte, rep_charge, rep_avg_mz, 
                                                             env_num, ec_score);
  SingleChargeFeaturePtrVec single_features;
  int charge = seed_ptr->getCharge();
  SingleChargeFeaturePtr single_feature = std::make_shared<SingleChargeFeature>(charge, ms1_time_begin, ms1_time_end,
                                                                                ms1_scan_begin,
                                                                                ms1_scan_end,
                                                                                feat_inte, env_num);
  single_features.push_back(single_feature);
  feature_ptr->setSingleFeatures(single_features);
  return feature_ptr;
}

}
}
