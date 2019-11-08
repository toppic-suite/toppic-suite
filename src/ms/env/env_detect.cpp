//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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
#include "common/base/mass_constant.hpp"
#include "ms/spec/raw_ms_util.hpp"
#include "ms/env/env_base.hpp"
#include "ms/env/env_detect.hpp"

namespace toppic {

MatchEnvPtr2D EnvDetect::getCandidate(const PeakPtrVec &peak_list, int max_charge, 
                                      double max_mass, EnvParaPtr env_para_ptr) {
  //PeakPtrVec peak_list = data_ptr->getPeakList();
  int peak_num = peak_list.size();
  //int max_charge = data_ptr->getMaxCharge();
  LOG_DEBUG("get candidate " << peak_num << " " << max_charge);
  MatchEnvPtr2D match_envs(peak_num);
  for (int idx = 0; idx < peak_num; idx++) {
    match_envs[idx].resize(max_charge);
    for (int charge = 1; charge <= max_charge; charge++) {
      match_envs[idx][charge - 1]
          = detectEnv(peak_list, idx, charge, max_mass, env_para_ptr);
    }
  }
  return match_envs;
}

// detect a MatchEnv
MatchEnvPtr EnvDetect::detectEnv(const PeakPtrVec &peak_list, int base_peak,
                                 int charge, double max_mass, EnvParaPtr env_para_ptr) {
  double base_mz = peak_list[base_peak]->getPosition();
  // check if the mass is greater than the precursor mass
  double base_mass = base_mz * charge - charge * mass_constant::getProtonMass();
  // LOG_DEBUG("peak idx " << base_peak << " base mz " << base_mz 
  // << " base mass " << base_mass << " max mass " << max_mass);
  if (base_mass >= max_mass || base_mass < env_para_ptr->min_mass_) {
    return nullptr;
  }

  // get a reference distribution based on the base mass
  EnvelopePtr ref_env = EnvBase::getStaticEnvByBaseMass(base_mass);
  if (ref_env == nullptr) {
    LOG_WARN("reference envelope is null");
    return nullptr;
  }

  // LOG_DEBUG("intensity " << peak_list[base_peak]->getIntensity() 
  // << " min refer inte " << env_para_ptr->min_refer_inte_);
  if (peak_list[base_peak]->getIntensity() < env_para_ptr->min_refer_inte_) {
    return nullptr;
  }
  // convert the reference distribution to a theoretical distribution
  // based on the base mz and charge state
  EnvelopePtr theo_env = ref_env->distrToTheoBase(base_mz, charge);
  int peak_idx = raw_ms_util::getNearPeakIdx(peak_list, theo_env->getReferMz(), 
                                             env_para_ptr->getMzTolerance());
  if (peak_idx < 0 || peak_list[peak_idx]->getIntensity() < env_para_ptr->min_refer_inte_) {
    return nullptr; 
  }
  // LOG_DEBUG("get theo env");
  // scale theoretical distribution
  double ratio = calcInteRatio(peak_list, theo_env, env_para_ptr->getMzTolerance());
  theo_env->changeIntensity(ratio);

  int mass_group = env_para_ptr->getMassGroup(base_mass);
  // LOG_DEBUG("theo env raw complete");

  // get the highest 85%--95% peaks
  double percentage = env_para_ptr->getPercentBound(mass_group);
  theo_env = theo_env->getSubEnv(percentage, env_para_ptr->min_inte_,
                                 env_para_ptr->max_back_peak_num_, env_para_ptr->max_forw_peak_num_);
  // get real envelope
  RealEnvPtr real_env = std::make_shared<RealEnv>(peak_list, theo_env, env_para_ptr->getMzTolerance(),
                                                  env_para_ptr->min_inte_);
  MatchEnvPtr match_env = std::make_shared<MatchEnv>(mass_group, theo_env, real_env);
  return match_env;
}

MatchEnvPtr EnvDetect::detectEnv(const PeakPtrVec &peak_list, double mono_mass,
                                 int charge, EnvParaPtr env_para_ptr) {
  if (mono_mass < env_para_ptr->min_mass_) {
    return nullptr;
  }

  // get a reference distribution based on the base mass
  EnvelopePtr ref_env = EnvBase::getStaticEnvByMonoMass(mono_mass);
  if (ref_env == nullptr) {
    LOG_WARN("reference envelope is null");
    return nullptr;
  }

  // convert the reference distribution to a theoretical distribution
  // based on the mono mz and charge state
  double mono_mz = mono_mass /charge + mass_constant::getProtonMass();
  EnvelopePtr theo_env = ref_env->distrToTheoMono(mono_mz, charge);
  // LOG_DEBUG("get theo env");
  // scale theoretical distribution
  double ratio = calcInteRatio(peak_list, theo_env, env_para_ptr->getMzTolerance());
  theo_env->changeIntensity(ratio);

  int mass_group = env_para_ptr->getMassGroup(mono_mass);
  // LOG_DEBUG("theo env raw complete");

  // get the highest 85%--95% peaks
  double percentage = env_para_ptr->getPercentBound(mass_group);
  theo_env = theo_env->getSubEnv(percentage, env_para_ptr->min_inte_,
                                 env_para_ptr->max_back_peak_num_, env_para_ptr->max_forw_peak_num_);
  // get real envelope
  RealEnvPtr real_env = std::make_shared<RealEnv>(peak_list, theo_env, env_para_ptr->getMzTolerance(),
                                                  env_para_ptr->min_inte_);
  MatchEnvPtr match_env = std::make_shared<MatchEnv>(mass_group, theo_env, real_env);
  return match_env;
}

// compute the intensity ratio based on the top three peaks
double EnvDetect::calcInteRatio(const PeakPtrVec &peak_list, EnvelopePtr theo_env,
                                double tolerance) {
  double theo_sum = 0;
  double obs_sum = 0;
  int refer_idx = theo_env->getReferIdx();
  double mz = theo_env->getMz(refer_idx);
  // LOG_DEBUG("step 1");
  int peak_idx = raw_ms_util::getNearPeakIdx(peak_list, mz, tolerance);
  if (peak_idx >= 0) {
    theo_sum += theo_env->getIntensity(refer_idx);
    obs_sum += peak_list[peak_idx]->getIntensity();
  }
  // LOG_DEBUG("step 2");
  if (refer_idx - 1 >= 0) {
    theo_sum += theo_env->getIntensity(refer_idx - 1);
    mz = theo_env->getMz(refer_idx-1);
    peak_idx = raw_ms_util::getNearPeakIdx(peak_list, mz, tolerance);
    if (peak_idx >= 0) {
      obs_sum += peak_list[peak_idx]->getIntensity();
    }
  }
  // LOG_DEBUG("step 3 refer idx " << refer_idx << " peak num " << theo_env->getPeakNum());
  if (refer_idx + 1 < theo_env->getPeakNum()) {
    theo_sum += theo_env->getIntensity(refer_idx + 1);
    mz = theo_env->getMz(refer_idx + 1);
    peak_idx = raw_ms_util::getNearPeakIdx(peak_list, mz, tolerance);
    if (peak_idx >= 0) {
      obs_sum += peak_list[peak_idx]->getIntensity();
    }
  }
  if (theo_sum == 0) {
    return 1.0;
  } else {
    return obs_sum / theo_sum;
  }
}

}  // namespace toppic
