// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "ms/env/match_env_refine.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <vector>

#include "common/base/mass_constant.hpp"
#include "common/util/logger.hpp"
#include "ms/env/env_base.hpp"

namespace toppic {

namespace match_env_refine {

void mzRefine(MatchEnvPtrVec& envs, double core_ratio) {
  for (size_t i = 0; i < envs.size(); i++) {
    mzRefine(envs[i], core_ratio);
  }
}

std::vector<int> coreIdxes(const ExpEnvPtr& real_env, const EnvPtr& theo_env,
                           double core_ratio) {
  std::vector<int> idxes;
  if (core_ratio <= 0) {
    return idxes;
  }
  double min_inte = core_ratio * theo_env->getReferInte();
  for (int i = 0; i < theo_env->getPeakNum(); i++) {
    if (theo_env->getInte(i) >= min_inte && real_env->isExist(i)) {
      idxes.push_back(i);
    }
  }
  if (idxes.size() < 3) {
    idxes.clear();
  }
  return idxes;
}

void mzRefine(const MatchEnvPtr& env, double core_ratio) {
  ExpEnvPtr real_env = env->getExpEnvPtr();
  double cur_mz = real_env->getReferMz();
  int charge = real_env->getCharge();
  double prev_mz = cur_mz - mass_constant::getIsotopeMass() / charge;
  double next_mz = cur_mz + mass_constant::getIsotopeMass() / charge;
  // check if the mass is greater than the precursor mass
  double ref_mass = peak_util::compPeakNeutralMass(cur_mz, charge);
  // get a reference distribution based on the reference mass
  EnvPtr refer_env = EnvBase::getEnvByRefMass(ref_mass);
  /* add one zeros at both sides of the envelope */
  EnvPtr ext_refer_env = refer_env->addZero(1);

  // convert the reference distribution to a theoretical distribution
  // based on the base mz and charge state
  int max_back_peak_num = real_env->getReferIdx();
  int max_forw_peak_num = real_env->getPeakNum() - real_env->getReferIdx() - 1;
  EnvPtr theo_env = ext_refer_env->distrToTheoRef(cur_mz, charge);
  double max_inte = theo_env->getReferInte();
  theo_env->changeIntensity(1.0 / max_inte);

  EnvPtr cur_env = theo_env->getSubEnv(max_back_peak_num, max_forw_peak_num);

  theo_env = ext_refer_env->distrToTheoRef(prev_mz, charge);
  max_inte = theo_env->getReferInte();
  theo_env->changeIntensity(1.0 / max_inte);
  EnvPtr prev_env;
  if (max_back_peak_num >= 1 &&
      real_env->isExist(real_env->getReferIdx() - 1)) {
    prev_env =
        theo_env->getSubEnv(max_back_peak_num - 1, max_forw_peak_num + 1);
  } else {
    prev_env = nullptr;
  }

  theo_env = ext_refer_env->distrToTheoRef(next_mz, charge);
  max_inte = theo_env->getReferInte();
  theo_env->changeIntensity(1.0 / max_inte);

  EnvPtr next_env;
  if (max_forw_peak_num >= 1 &&
      real_env->isExist(real_env->getReferIdx() + 1)) {
    next_env =
        theo_env->getSubEnv(max_back_peak_num + 1, max_forw_peak_num - 1);
  } else {
    next_env = nullptr;
  }
  // Stage 1: choose the monoisotopic position (current, previous or next)
  // with the whole-envelope distance, as before. On a myoglobin MS/MS
  // spectrum a core-only choice matched no more theoretical fragment masses
  // than this one, so the choice is left unchanged.
  std::vector<int> all_idxes;  // empty = all peaks
  double cur_dist;
  double cur_ratio;
  compEnvDist(real_env, cur_env, all_idxes, cur_dist, cur_ratio);
  double prev_dist;
  double prev_ratio;
  compEnvDist(real_env, prev_env, all_idxes, prev_dist, prev_ratio);
  double next_dist;
  double next_ratio;
  compEnvDist(real_env, next_env, all_idxes, next_dist, next_ratio);

  EnvPtr chosen_env;
  if (cur_dist <= prev_dist && cur_dist <= next_dist) {
    chosen_env = cur_env;
  } else if (prev_dist <= next_dist) {
    int peak_num = prev_env->getPeakNum();
    if (prev_env->getInte(peak_num - 1) == 0) {
      prev_env->removeRightPeaks(1);
      real_env->removeRightPeaks(1);
    }
    chosen_env = prev_env;
    real_env->changeReferIdx(-1);
  } else {
    if (next_env->getInte(0) == 0) {
      next_env->removeLeftPeaks(1);
      real_env->removeLeftPeaks(1);
    }
    chosen_env = next_env;
    real_env->changeReferIdx(1);
  }

  // Stage 2: scale the chosen distribution by refitting its intensity ratio
  // on the core peaks only. The whole-envelope ratio is pulled up when the
  // envelope's tails are inflated by overlapping neighbouring envelopes,
  // leaving a theoretical apex well above the observed one; the core (peaks
  // near the reference peak) is the part least affected by such overlap.
  // coreIdxes returns an empty list (= all peaks, i.e. the stage-1 ratio) for
  // envelopes with fewer than 3 usable core peaks.
  std::vector<int> core_idxes = coreIdxes(real_env, chosen_env, core_ratio);
  double core_dist;
  double chosen_ratio;
  compEnvDist(real_env, chosen_env, core_idxes, core_dist, chosen_ratio);
  chosen_env->changeIntensity(chosen_ratio);
  env->setTheoEnvPtr(chosen_env);
}

void compEnvDist(const EnvPtr& real_env, const EnvPtr& theo_env,
                 const std::vector<int>& idxes, double& dist, double& ratio) {
  if (theo_env == nullptr) {
    dist = std::numeric_limits<double>::infinity();
  } else {
    compDistWithNorm(real_env->getInteList(), theo_env->getInteList(), idxes,
                     dist, ratio);
  }
}

void compDistWithNorm(const std::vector<double>& real,
                      const std::vector<double>& theo,
                      const std::vector<int>& idxes, double& best_dist,
                      double& best_ratio) {
  best_dist = std::numeric_limits<double>::infinity();
  best_ratio = -1;
  // Empty idxes means every peak.
  std::vector<int> all_idxes;
  const std::vector<int>* fit_idxes = &idxes;
  if (idxes.empty()) {
    all_idxes.resize(real.size());
    std::iota(all_idxes.begin(), all_idxes.end(), 0);
    fit_idxes = &all_idxes;
  }
  // Candidate ratios come from the fitted peaks only, so a contaminated tail
  // peak can neither seed nor score the search.
  for (int i : *fit_idxes) {
    if (theo[i] == 0.0) {
      continue;
    }
    double ratio = real[i] / theo[i];
    if (ratio <= 0) {
      continue;
    }
    for (int j = 80; j <= 120; j++) {
      double cur_ratio = ratio * j / 100;
      std::vector<double> norm_real = norm(real, cur_ratio);
      double dist = compDist(norm_real, theo, *fit_idxes);
      if (dist < best_dist) {
        best_dist = dist;
        best_ratio = cur_ratio;
      }
    }
  }
  if (std::isnan(best_ratio)) {
    LOG_ERROR("The best ratio is not a number!");
  }
}

std::vector<double> norm(const std::vector<double>& obs, double ratio) {
  std::vector<double> result(obs.size());
  for (size_t i = 0; i < obs.size(); i++) {
    result[i] = obs[i] / ratio;
  }
  return result;
}

double compDist(const std::vector<double>& norm,
                const std::vector<double>& theo,
                const std::vector<int>& idxes) {
  double max_distance_a = 1.0;
  double max_distance_b = 1.0;
  double result = 0;
  std::vector<int> all_idxes;
  const std::vector<int>* fit_idxes = &idxes;
  if (idxes.empty()) {
    all_idxes.resize(norm.size());
    std::iota(all_idxes.begin(), all_idxes.end(), 0);
    fit_idxes = &all_idxes;
  }
  for (int i : *fit_idxes) {
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
