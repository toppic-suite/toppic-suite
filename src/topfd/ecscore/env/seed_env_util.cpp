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

#include <cmath>

#include "common/util/logger.hpp"
#include "common/base/mass_constant.hpp"
#include "topfd/ecscore/env/ms_map_env_util.hpp"
#include "topfd/ecscore/env/seed_env_util.hpp"

namespace toppic {

namespace seed_env_util {

bool containEnoughPeaks(EnvPeakPtrVec &peak_ptr_list, int charge) {
  int peak_num = 0;
  for (size_t i = 0; i < peak_ptr_list.size(); i++) {
    if (peak_ptr_list[i]->isExist()) {
      peak_num++;
    }
  }
  if ((charge == 1 || charge == 2) && peak_num < 2) return false;
  if (charge > 2 && charge < 15 && peak_num < 3) return false;
  if (charge >= 15 && peak_num < 5) return false;
  return true;
}

double pearsonr(std::vector<double> &X, std::vector<double> &Y) {
  int n = X.size();
  double sum_X = 0, sum_Y = 0, sum_XY = 0;
  double squareSum_X = 0, squareSum_Y = 0;
  for (int i = 0; i < n; i++) {
    sum_X = sum_X + X[i];
    sum_Y = sum_Y + Y[i];
    sum_XY = sum_XY + X[i] * Y[i];
    squareSum_X = squareSum_X + X[i] * X[i];
    squareSum_Y = squareSum_Y + Y[i] * Y[i];
  }
  double corr = (double) (n * sum_XY - sum_X * sum_Y) /
    sqrt((n * squareSum_X - sum_X * sum_X) * (n * squareSum_Y - sum_Y * sum_Y));
  return corr;
}

double pearsonr(EnvPeakPtrVec &peak_list, std::vector<double> &Y) {
  int n = 0;
  double sum_X = 0, sum_Y = 0, sum_XY = 0;
  double squareSum_X = 0, squareSum_Y = 0;
  for (size_t i = 0; i < peak_list.size(); i++) {
    if (!peak_list[i]->isExist()) {
      continue;
    }
    double x = peak_list[i]->getIntensity();
    sum_X = sum_X + x;
    sum_Y = sum_Y + Y[i];
    sum_XY = sum_XY + x * Y[i];
    squareSum_X = squareSum_X + x * x;
    squareSum_Y = squareSum_Y + Y[i] * Y[i];
    n++;
  }
  double corr = (double) (n * sum_XY - sum_X * sum_Y) /
    sqrt((n * squareSum_X - sum_X * sum_X) * (n * squareSum_Y - sum_Y * sum_Y));
  return corr;
}


SeedEnvPtr preprocessSeedEnvPtr(SeedEnvPtr seed_ptr, MsMapPtr ms_map_ptr,
                                EcscoreParaPtr para_ptr, double sn_ratio) {
  //1. Check charge state
  if (seed_ptr->getCharge() < para_ptr->para_min_charge_) {
    return nullptr;
  }
  //2. Check spectrum id
  if (seed_ptr->getSpecId() >= ms_map_ptr->getRowNum()) {
    LOG_ERROR("spec id " + std::to_string(seed_ptr->getSpecId()) + " is out of range!");
    return nullptr;
  }
  //3. Get high intensity peaks only. 
  MsMapEnvPtr ms_map_env_ptr
    = ms_map_env_util::getMatchMsMapEnv(ms_map_ptr, seed_ptr,
                                        seed_ptr->getSpecId(),
                                        para_ptr->mass_tole_);
  double inte_ratio = ms_map_env_util::compTopThreeInteRatio(seed_ptr, ms_map_env_ptr);
  double min_inte = ms_map_ptr->getBaseInte() * sn_ratio;
  EnvPeakPtrVec scaled_peak_ptr_list = seed_ptr->getScaledPeakPtrList(inte_ratio, min_inte);
  //4. Check peak number
  if (!containEnoughPeaks(scaled_peak_ptr_list, seed_ptr->getCharge())){
    return nullptr;
  }
  //5. Check Pearson correlation 
  std::vector<double> exp_inte_list = ms_map_env_ptr->getInteList();
  double corr = pearsonr(scaled_peak_ptr_list, exp_inte_list); 
  if (corr < para_ptr->corr_tole_) {
    return nullptr;
  }
  //6. Generate new seed envelope
  SeedEnvPtr result_seed_env_ptr = std::make_shared<SeedEnv>(seed_ptr,
                                                             scaled_peak_ptr_list); 
  return result_seed_env_ptr;
}

bool containOnePeak(EnvPeakPtrVec &peak_ptr_list) {
  int peak_num = 0;
  for (size_t i = 0; i < peak_ptr_list.size(); i++) {
    if (peak_ptr_list[i]->isExist()) {
      peak_num++;
    }
  }
  if (peak_num >= 1) {
    return true;
  }
  else {
    return false;
  }
}

SeedEnvPtr relaxProcessSeedEnvPtr(SeedEnvPtr seed_ptr, MsMapPtr ms_map_ptr,
                                  EcscoreParaPtr para_ptr, double sn_ratio) {
  //1. Check charge state
  if (seed_ptr->getCharge() < para_ptr->para_min_charge_) {
    return nullptr;
  }
  //2. Check spectrum id
  if (seed_ptr->getSpecId() >= ms_map_ptr->getRowNum()) {
    LOG_ERROR("spec id " + std::to_string(seed_ptr->getSpecId()) + " is out of range!");
    return nullptr;
  }
  //3. Get high intensity peaks only. 
  MsMapEnvPtr ms_map_env_ptr
    = ms_map_env_util::getMatchMsMapEnv(ms_map_ptr, seed_ptr,
                                        seed_ptr->getSpecId(),
                                        para_ptr->mass_tole_);
  double inte_ratio = ms_map_env_util::compTopThreeInteRatio(seed_ptr, ms_map_env_ptr);
  double min_inte = ms_map_ptr->getBaseInte() * sn_ratio;
  EnvPeakPtrVec scaled_peak_ptr_list = seed_ptr->getScaledPeakPtrList(inte_ratio, min_inte);
  //4. Check peak number
  if (!containOnePeak(scaled_peak_ptr_list)){
    return nullptr;
  }
  //5. Generate new seed envelope
  SeedEnvPtr result_seed_env_ptr = std::make_shared<SeedEnv>(seed_ptr,
                                                             scaled_peak_ptr_list); 
  return result_seed_env_ptr;
}

SeedEnvPtr getHalfChargeEnvV1(SeedEnvPtr seed_ptr,
                              double even_odd_peak_ratio) {
  double mass = seed_ptr->getMonoNeutralMass();
  int charge = seed_ptr->getCharge();
  double mz = peak_util::compMz(mass, charge);
  std::vector<double> distribution = seed_ptr->getMzList();
  if (even_odd_peak_ratio < 0) {
    mz = mz + (distribution[1] - distribution[0]);
  }
  int new_charge = int(charge / 2);
  if (new_charge == 0) {
    new_charge = new_charge + 1;
  }
  double mono_mass = peak_util::compPeakNeutralMass(mz, new_charge);

  int sp_id = seed_ptr->getSpecId();
  int peak_id = -1;
  double inte = seed_ptr->getSeedInte()/2;

  DeconvPeakPtr peak_ptr = std::make_shared<DeconvPeak>(sp_id, peak_id,
                                                        mono_mass, inte,
                                                        new_charge);
  SeedEnvPtr new_seed_ptr = std::make_shared<SeedEnv>(peak_ptr);

  return new_seed_ptr;
}

SeedEnvPtr getHalfChargeEnvV2(SeedEnvPtr seed_ptr,
                              double even_odd_log_ratio) {
  double old_charge = seed_ptr->getCharge();
  if (old_charge < 2) {
    return nullptr;
  }
  int new_charge = int(old_charge / 2);
  double ref_mz = seed_ptr->getReferMz();
  double ref_mass = peak_util::compPeakNeutralMass(ref_mz, new_charge);
  double mono_mass = EnvBase::convertRefMassToMonoMass(ref_mass);
  int refer_idx = seed_ptr->getReferIdx();
  // if refer_idx is even and env_odd_log_ratio < 0 
  // or refer_idx is odd  and env_odd_log_ratio > 1
  // then increase mass by about 1 Dalton
  if (((refer_idx)%2 == 0 && even_odd_log_ratio < 0) 
      || ((refer_idx)%2 == 1 && even_odd_log_ratio >0)) {
    mono_mass += mass_constant::getIsotopeMass();
  }
  int sp_id = seed_ptr->getSpecId();
  int peak_id = -1;
  double inte = seed_ptr->getSeedInte()/2;
  DeconvPeakPtr peak_ptr = std::make_shared<DeconvPeak>(sp_id, peak_id,
                                                        mono_mass, inte,
                                                        new_charge);
  SeedEnvPtr new_seed_ptr = std::make_shared<SeedEnv>(peak_ptr);
  return new_seed_ptr;
}

SeedEnvPtr testHalfChargeEnv(SeedEnvPtr seed_ptr, MsMapPtr ms_map_ptr, 
                             double even_odd_log_ratio, EcscoreParaPtr para_ptr, 
                             double sn_ratio) {
  SeedEnvPtr half_charge_seed = getHalfChargeEnvV1(seed_ptr, even_odd_log_ratio);
  if (half_charge_seed == nullptr) {
      return nullptr;
  }
  SeedEnvPtr processed_seed_ptr = preprocessSeedEnvPtr(half_charge_seed, 
                                                       ms_map_ptr, para_ptr, sn_ratio);
  return processed_seed_ptr;
}


}

}
