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

#include <limits>
#include <cmath>

#include "common/util/logger.hpp"
#include "common/base/mass_constant.hpp"
#include "deconv/env/env_base.hpp"
#include "deconv/env/match_env_refine.hpp"

namespace toppic {

namespace match_env_refine {

void mzRefine(MatchEnvPtrVec &envs) {
  for (size_t i = 0; i < envs.size(); i++) {
    mzRefine(envs[i]);
  }
}

void mzRefine(MatchEnvPtr env) {
  RealEnvPtr real_env = env->getRealEnvPtr();
  double cur_mz = real_env->getReferMz();
  int charge = real_env->getCharge();
  double prev_mz = cur_mz - 1.0 / charge;
  double next_mz = cur_mz + 1.0 / charge;
  // check if the mass is greater than the precursor mass
  double bass_mass = cur_mz * charge - charge * mass_constant::getProtonMass();
  // get a reference distribution based on the base mass
  EnvelopePtr refer_env = EnvBase::getStaticEnvByBaseMass(bass_mass);
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
  double cur_dist;
  double cur_ratio;
  compEnvDist(real_env, cur_env, cur_dist, cur_ratio);
  double prev_dist;
  double prev_ratio;
  compEnvDist(real_env, prev_env, prev_dist, prev_ratio);
  double next_dist;
  double next_ratio;
  compEnvDist(real_env, next_env, next_dist, next_ratio);

  if (cur_dist <= prev_dist && cur_dist <= next_dist) {
    cur_env->changeIntensity(cur_ratio);
    env->setTheoEnvPtr(cur_env);
  } else if (prev_dist <= next_dist) {
    prev_env->changeIntensity(prev_ratio);
    env->setTheoEnvPtr(prev_env);
    real_env->shift(-1);
  } else {
    next_env->changeIntensity(next_ratio);
    env->setTheoEnvPtr(next_env);
    real_env->shift(1);
  }
}

void compEnvDist(EnvelopePtr real_env, EnvelopePtr theo_env, 
                 double &dist, double &ratio) {
  if (theo_env == nullptr) {
    dist = std::numeric_limits<double>::infinity();
  } else {
    compDistWithNorm(real_env->getIntensities(), theo_env->getIntensities(),
                     dist, ratio);
  }
}

void compDistWithNorm(const std::vector<double>& real,
                      const std::vector<double>& theo, 
                      double &best_dist,
                      double &best_ratio) {
  best_dist = std::numeric_limits<double>::infinity();
  best_ratio = -1;
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
        best_ratio = cur_ratio;
      }
    }
  }
}

std::vector<double> norm(const std::vector<double> &obs, double ratio) {
  std::vector<double> result(obs.size());
  for (size_t i = 0; i < obs.size(); i++) {
    result[i] = obs[i] / ratio;
  }
  return result;
}

double compDist(const std::vector<double> &norm, const std::vector<double> &theo) {
  double max_distance_a = 1.0;
  double max_distance_b = 1.0;
  double result = 0;
  for (size_t i = 0; i < norm.size(); i++) {
    double dist = std::abs(norm[i] - theo[i]);
    if (norm[i] > theo[i]) {
      if (dist > max_distance_a) {
        dist = max_distance_a;
      }
    } else {
      if (dist > max_distance_b) {
        dist = max_distance_b;
      }
    }
    result = result + dist * dist;
  }
  return result;
}

}  // namespace match_env_refine

}  // namespace toppic
