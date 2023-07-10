//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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
#include "ms/spec/peak_list_util.hpp"
#include "ms/env/env_base.hpp"
#include "ms/env/env_detect.hpp"

namespace toppic {

namespace env_detect {

// compute the intensity ratio based on the top three peaks
double calcInteRatio(const PeakPtrVec &peak_list, EnvelopePtr theo_env_ptr,
                     double tolerance) {
  double theo_sum = 0;
  double obs_sum = 0;
  int refer_idx = theo_env_ptr->getReferIdx();
  double mz = theo_env_ptr->getMz(refer_idx);
  int peak_idx = peak_list_util::getNearPeakIdx(peak_list, mz, tolerance);
  if (peak_idx >= 0) {
    theo_sum += theo_env_ptr->getIntensity(refer_idx);
    obs_sum += peak_list[peak_idx]->getIntensity();
  }
  if (refer_idx - 1 >= 0) {
    theo_sum += theo_env_ptr->getIntensity(refer_idx - 1);
    mz = theo_env_ptr->getMz(refer_idx-1);
    peak_idx = peak_list_util::getNearPeakIdx(peak_list, mz, tolerance);
    if (peak_idx >= 0) {
      obs_sum += peak_list[peak_idx]->getIntensity();
    }
  }
  if (refer_idx + 1 < theo_env_ptr->getPeakNum()) {
    theo_sum += theo_env_ptr->getIntensity(refer_idx + 1);
    mz = theo_env_ptr->getMz(refer_idx + 1);
    peak_idx = peak_list_util::getNearPeakIdx(peak_list, mz, tolerance);
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

MatchEnvPtr compEnv(const PeakPtrVec &peak_list, EnvelopePtr theo_env_ptr, 
                    EnvParaPtr env_para_ptr, int mass_group, int min_inte) {
  // scale theoretical distribution
  double ratio = calcInteRatio(peak_list, theo_env_ptr, env_para_ptr->getMzTolerance());
  theo_env_ptr->changeIntensity(ratio);

  // get the highest 85%--95% peaks
  double percentage = env_para_ptr->getPercentBound(mass_group);
  theo_env_ptr = theo_env_ptr->getSubEnv(percentage, min_inte,
                                         env_para_ptr->max_back_peak_num_, 
                                         env_para_ptr->max_forw_peak_num_);
  // get real envelope
  RealEnvPtr real_env_ptr = std::make_shared<ExpEnv>(peak_list, theo_env_ptr,
                                                     env_para_ptr->getMzTolerance(),
                                                     min_inte);
  MatchEnvPtr match_env_ptr = std::make_shared<MatchEnv>(mass_group, theo_env_ptr, real_env_ptr);
  return match_env_ptr;
}

// detect a MatchEnv
MatchEnvPtr detectEnvByRefPeak(const PeakPtrVec &peak_list, int ref_peak, int charge, 
                               double max_mass, double min_inte, double min_ref_inte, 
                               EnvParaPtr env_para_ptr) {

  double refer_mz = peak_list[ref_peak]->getPosition();
  // check if the mass is greater than the precursor mass
  double ref_mass = peak_util::compPeakNeutralMass(refer_mz, charge); 
  if (ref_mass >= max_mass || ref_mass < env_para_ptr->min_mass_) {
    return nullptr;
  }

  // get a reference distribution based on the reference (highest intensity) mass
  EnvelopePtr ref_env_ptr = EnvBase::getEnvByRefMass(ref_mass);
  if (ref_env_ptr == nullptr) {
    LOG_WARN("reference envelope is null");
    return nullptr;
  }

  if (peak_list[ref_peak]->getIntensity() < min_ref_inte) {
    return nullptr;
  }
  // convert the distribution to a theoretical distribution
  // based on the refer mz and charge state
  EnvelopePtr theo_env_ptr = ref_env_ptr->distrToTheoRef(refer_mz, charge);
  int peak_idx = peak_list_util::getNearPeakIdx(peak_list, theo_env_ptr->getReferMz(), 
                                                env_para_ptr->getMzTolerance());
  if (peak_idx < 0 || peak_list[peak_idx]->getIntensity() < min_ref_inte) {
    return nullptr; 
  }
  int mass_group = env_para_ptr->getMassGroup(ref_mass);

  return compEnv(peak_list, theo_env_ptr, env_para_ptr, mass_group, min_inte);
}

MatchEnvPtr detectEnvByMonoMass(const PeakPtrVec &peak_list, double mono_mass,
                                int charge, double min_inte, EnvParaPtr env_para_ptr) {
  if (mono_mass < env_para_ptr->min_mass_) {
    return nullptr;
  }

  // get a reference distribution based on the base mass
  EnvelopePtr ref_env_ptr = EnvBase::getEnvByMonoMass(mono_mass);
  if (ref_env_ptr == nullptr) {
    LOG_WARN("Reference envelope is null!");
    return nullptr;
  }

  // convert the reference distribution to a theoretical distribution
  // based on the mono mz and charge state
  double mono_mz = mono_mass /charge + mass_constant::getProtonMass();
  EnvelopePtr theo_env_ptr = ref_env_ptr->distrToTheoMono(mono_mz, charge);

  int mass_group = env_para_ptr->getMassGroup(mono_mass);

  return compEnv(peak_list, theo_env_ptr, env_para_ptr, mass_group, min_inte);
}

MatchEnvPtr2D getCandidateEnv(const PeakPtrVec &peak_list, int max_charge, double max_mass, 
                              double min_inte, double min_ref_inte, EnvParaPtr env_para_ptr) {
  int peak_num = peak_list.size();
  MatchEnvPtr2D match_envs(peak_num);
  for (int idx = 0; idx < peak_num; idx++) {
    match_envs[idx].resize(max_charge);
    for (int charge = 1; charge <= max_charge; charge++) {
      match_envs[idx][charge - 1]
          = detectEnvByRefPeak(peak_list, idx, charge, max_mass, min_inte, min_ref_inte, env_para_ptr);
    }
  }
  return match_envs;
}

}  // namespace env_detect

}  // namespace toppic
