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


#include "common/util/logger.hpp"
#include "topfd/ecscore/env/env_util.hpp"
#include "topfd/ecscore/env/seed_env_util.hpp"
#include "topfd/ecscore/env_set/env_set_util.hpp"

namespace toppic {

namespace seed_env_util {

bool preprocessEnv(MsMapPtr matrix_ptr, SeedEnvPtr seed_ptr,
                   EcscoreParaPtr para_ptr, double sn_ratio) {
  if (seed_ptr->getCharge() < para_ptr->para_min_charge_)
    return false;
  double mass_tol = para_ptr->mass_tole_;
  double corr_tol = para_ptr->corr_tole_;
  //double min_mz = matrix_ptr->getMinMz() - mass_tol;
  //double max_mz = matrix_ptr->getMaxMz() + mass_tol;
  //seed_ptr->seedRmPeaks(min_mz, max_mz);
  //env_set_util::compPeakStartEndIdx(matrix_ptr, seed_ptr, mass_tol);
  bool valid = evalEnv(matrix_ptr, seed_ptr, mass_tol, corr_tol, sn_ratio);
  if (seed_ptr->getSpecId() >= matrix_ptr->getRowNum()) {
    LOG_ERROR("spec id " + std::to_string(seed_ptr->getSpecId()) + " is out of range!");
    valid = false;
  }
  return valid;
}

bool evalEnv(MsMapPtr matrix_ptr, SeedEnvPtr seed_ptr,
             double mass_tol, double corr_tol, double sn_ratio){
  double noise_inte = matrix_ptr->getBaseInte();
  ExpEnvelopePtr exp_env_ptr 
    = env_set_util::getMatchExpEnv(matrix_ptr, seed_ptr, seed_ptr->getSpecId(), mass_tol);
  std::vector<double> exp_env_mass = exp_env_ptr->getPosList();
  std::vector<double> seed_env_mass = seed_ptr->getMzList();
  std::vector<double> exp_env_inte = exp_env_ptr->getInteList();
  std::vector<double> seed_env_inte = seed_ptr->getInteList();
  int num_peaks = seed_env_inte.size();
  double inte_ratio = env_set_util::calcInteRatio(seed_env_inte, exp_env_inte);
  std::vector<double> scaled_theo_inte;
  for (int i = 0; i < num_peaks; i++) {
    double scaled_inte = inte_ratio * seed_env_inte[i];
    if (scaled_inte < sn_ratio * noise_inte)
      scaled_inte = 0;
    scaled_theo_inte.push_back(scaled_inte);
  }
  EnvPeakPtrVec seed_env_peaks = seed_ptr->getPeakPtrList();
  for (int j = num_peaks-1; j >= 0; j--) {
    if (scaled_theo_inte[j] == 0) {
      scaled_theo_inte.erase(scaled_theo_inte.begin() + j);
      seed_env_peaks.erase(seed_env_peaks.begin() + j);
      exp_env_inte.erase(exp_env_inte.begin() + j);
    }
  }
  seed_ptr->setPeakPtrList(seed_env_peaks);
  if (!testChargeState(seed_ptr->getCharge(), scaled_theo_inte)) {
    return false;
  }
  return evalEnvPair(exp_env_inte, scaled_theo_inte, corr_tol);
}

bool simpleEvalEnv(MsMapPtr matrix_ptr, SeedEnvPtr seed_ptr,
                   double mass_tol, double corr_tol, double sn_ratio){
  int spec_id = seed_ptr->getSpecId();
  ExpEnvelopePtr exp_env_ptr 
    = env_set_util::getMatchExpEnv(matrix_ptr, seed_ptr, spec_id, mass_tol);
  std::vector<double> exp_env_mass = exp_env_ptr->getPosList();
  std::vector<double> seed_env_mass = seed_ptr->getMzList();
  std::vector<double> exp_env_inte = exp_env_ptr->getInteList();
  std::vector<double> seed_env_inte = seed_ptr->getInteList();
  int num_peaks = seed_env_inte.size();
  double inte_ratio = env_set_util::calcInteRatio(seed_env_inte, exp_env_inte);
  double noise_inte = matrix_ptr->getBaseInte();
  std::vector<double> scaled_theo_inte;
  for (int i = 0; i < num_peaks; i++) {
    double scaled_inte = inte_ratio * seed_env_inte[i];
    if (scaled_inte < sn_ratio * noise_inte)
      scaled_inte = 0;
    scaled_theo_inte.push_back(scaled_inte);
  }
  EnvPeakPtrVec seed_env_peaks = seed_ptr->getPeakPtrList();
  for (int j = num_peaks-1; j >= 0; j--) {
    if (scaled_theo_inte[j] == 0) {
      scaled_theo_inte.erase(scaled_theo_inte.begin() + j);
      seed_env_peaks.erase(seed_env_peaks.begin() + j);
      exp_env_inte.erase(exp_env_inte.begin() + j);
    }
  }
  seed_ptr->setPeakPtrList(seed_env_peaks);
  if (seed_env_inte.size() < 1) {
    return false;
  }
  int num_theo_peak = 0;
  for (auto i : scaled_theo_inte) { 
    if (i > 0) {
      num_theo_peak++;
    }
  }
  if (num_theo_peak == 0) return false;
  return true;
}

bool simplePreprocessEnv(MsMapPtr matrix_ptr, SeedEnvPtr seed_ptr,
                         EcscoreParaPtr para_ptr, double sn_ratio) {
  if (seed_ptr->getSpecId() >= matrix_ptr->getRowNum()) {
    LOG_ERROR("spec id " + std::to_string(seed_ptr->getSpecId()) + " is out of range!");
    return false; 
  }
  if (seed_ptr->getCharge() < para_ptr->para_min_charge_) {
    return false;
  }
  double mass_tol = para_ptr->mass_tole_;
  double corr_tol = para_ptr->corr_tole_;
  //double min_mz = matrix_ptr->getMinMz() - mass_tol;
  //double max_mz = matrix_ptr->getMaxMz() + mass_tol;
  std::vector<double> seed_env_inte = seed_ptr->getInteList();
  //  seed_ptr->seedRmPeaks(min_mz, max_mz);
  seed_env_inte = seed_ptr->getInteList();
  //env_set_util::compPeakStartEndIdx(matrix_ptr, seed_ptr, mass_tol);
  seed_env_inte = seed_ptr->getInteList();
  bool valid = simpleEvalEnv(matrix_ptr, seed_ptr, mass_tol, corr_tol, sn_ratio);
  return valid;
}


bool testChargeState(int charge, std::vector<double> &seed_env_inte) {
  if ((charge == 1 || charge == 2) && seed_env_inte.size() < 2) return false;
  if (charge > 2 && charge < 15 && seed_env_inte.size() < 3) return false;
  if (charge >= 15 && seed_env_inte.size() < 5) return false;
  return true;
}


bool evalEnvPair(std::vector<double> &exp_env_inte, 
                 std::vector<double> &theo_inte, double tol) {
  int num_theo_peak = 0;
  for (auto i : theo_inte) { 
    if (i > 0) {
      num_theo_peak++;
    }
  }
  if (num_theo_peak == 0) return false;
  double corr = env_util::pearsonr(exp_env_inte, theo_inte);
  if (corr < tol) { 
    return false;
  }
  else {
    return true;
  }
}

}

}
