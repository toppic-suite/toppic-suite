// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include <algorithm>

#include "base/logger.hpp"
#include "feature/match_env_refine.hpp" 

namespace prot {

double MatchEnvRefine::max_distance_a_ = 1.0;
double MatchEnvRefine::max_distance_b_ = 1.0;
double MatchEnvRefine::best_ratio_;

void MatchEnvRefine::mzRefine(FeatureMngPtr mng_ptr, MatchEnvPtrVec &envs) {
  for (size_t i = 0; i < envs.size(); i++) {
    mzRefine(mng_ptr, envs[i]);
  }
}

void MatchEnvRefine::mzRefine(FeatureMngPtr mng_ptr, MatchEnvPtr env) {
  RealEnvPtr real_env = env->getRealEnvPtr();
  double cur_mz = real_env->getReferMz();
  int charge = real_env->getCharge();
  double prev_mz = cur_mz - 1.0 / charge;
  double next_mz = cur_mz + 1.0 / charge;
  // check if the mass is greater than the precursor mass 
  double bass_mass = cur_mz * charge - charge * MassConstant::getProtonMass();
  // get a reference distribution based on the base mass 
  EnvelopePtr refer_env = mng_ptr->env_base_ptr_->getEnvByBaseMass(bass_mass);
  /* add one zeros at both sides of the envelope */
  EnvelopePtr ext_refer_env = refer_env->addZero(1);

  // convert the reference distribution to a theoretical distribution
  // based on the base mz and charge state
  int max_back_peak_num = real_env->getReferIdx();
  int max_forw_peak_num = real_env->getPeakNum() - real_env->getReferIdx() - 1;
  EnvelopePtr theo_env = ext_refer_env->distrToTheoBase(cur_mz, charge);
  double max_inte = theo_env->getReferIntensity();
  theo_env->changeIntensity(1.0 / max_inte);

  EnvelopePtr cur_env = theo_env->getSubEnv(max_back_peak_num, max_forw_peak_num);

  theo_env = ext_refer_env->distrToTheoBase(prev_mz, charge);
  max_inte = theo_env->getReferIntensity();
  theo_env->changeIntensity(1.0 / max_inte);
  EnvelopePtr prev_env;
  if (max_back_peak_num >= 1 && real_env->isExist(real_env->getReferIdx() - 1)) {
    prev_env = theo_env->getSubEnv(max_back_peak_num - 1, max_forw_peak_num + 1);
  } else {
    prev_env = nullptr;
  }

  theo_env = ext_refer_env->distrToTheoBase(next_mz, charge);
  max_inte = theo_env->getReferIntensity();
  theo_env->changeIntensity(1.0 / max_inte);

  EnvelopePtr next_env;
  if (max_forw_peak_num >= 1 && real_env->isExist(real_env->getReferIdx() + 1)) {
    next_env = theo_env->getSubEnv(max_back_peak_num + 1, max_forw_peak_num - 1);
  } else {
    next_env = nullptr;
  }

  double cur_dist = compEnvDist(real_env, cur_env);
  double cur_ratio = best_ratio_;
  double prev_dist = compEnvDist(real_env, prev_env);
  double prev_ratio = best_ratio_;
  double next_dist = compEnvDist(real_env, next_env);
  double next_ratio = best_ratio_;

  if (cur_dist <= prev_dist && cur_dist <= next_dist) {
    cur_env->changeIntensity(cur_ratio);
    env->setTheoEnvPtr(cur_env);
  } else if (prev_dist <= next_dist){
    prev_env->changeIntensity(prev_ratio);
    env->setTheoEnvPtr(prev_env);
    real_env->shift(-1);
  }
  else {
    next_env->changeIntensity(next_ratio);
    env->setTheoEnvPtr(next_env);
    real_env->shift(1);
  }
}

double MatchEnvRefine::compEnvDist(EnvelopePtr real_env, EnvelopePtr theo_env) {
  if (theo_env == nullptr) {
    return std::numeric_limits<double>::infinity();
  } else {
    return compDistWithNorm(real_env->getIntensities(), theo_env->getIntensities());
  }
}

double MatchEnvRefine::compDistWithNorm(std::vector<double> real, std::vector<double> theo) {
  double best_dist = std::numeric_limits<double>::infinity();
  best_ratio_ = -1;
  for (size_t i = 0; i < real.size(); i++) {
    double ratio = real[i] / theo[i];
    if (ratio <= 0) {
      continue;
    }
    for (int j = 80; j <= 120; j++) {
      double cur_ratio = ratio * j / 100;
      std::vector<double> norm_real = norm(real, cur_ratio);
      double dist = compDist(norm_real, theo);
      if (dist < best_dist) {
        best_dist = dist;
        best_ratio_ = cur_ratio;
      }
    }
  }
  return best_dist;
}

std::vector<double> MatchEnvRefine::norm(std::vector<double> &obs, double ratio) {
  std::vector<double> result;
  for (size_t i = 0; i < obs.size(); i++) {
    result.push_back(obs[i] / ratio);
  }
  return result;
}

double MatchEnvRefine::compDist(std::vector<double> &norm, std::vector<double> &theo) {
  double result = 0;
  for (size_t i = 0; i < norm.size(); i++) {
    double dist = std::abs(norm[i] - theo[i]);
    if (norm[i] > theo[i]) {
      if (dist > max_distance_a_) {
        dist = max_distance_a_;
      }
    } else {
      if (dist > max_distance_b_) {
        dist = max_distance_b_;
      }
    }
    result = result + dist * dist;
  }
  return result;
}

}
