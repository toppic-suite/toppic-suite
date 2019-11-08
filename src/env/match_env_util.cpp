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


#include <algorithm>

#include "common/util/logger.hpp"
#include "common/base/mass_constant.hpp"
#include "spec/raw_ms_util.hpp"
#include "env/match_env_util.hpp"

namespace toppic {

namespace match_env_util {

std::vector<double> getMassList(const MatchEnvPtrVec &envs) {
  std::vector<double> masses(envs.size());
  for (size_t i = 0; i < envs.size(); i++) {
    masses[i] = envs[i]->getRealEnvPtr()->getMonoNeutralMass();
  }
  return masses;
}

std::vector<int> getChargeList(const MatchEnvPtrVec &envs) {
  std::vector<int> charges(envs.size());
  for (size_t i = 0; i < envs.size(); i++) {
    charges[i] = envs[i]->getRealEnvPtr()->getCharge();
  }
  return charges;
}

std::vector<double> getChargeOneMassList(const MatchEnvPtrVec &envs) {
  std::vector<double> masses(envs.size());
  for (size_t i = 0; i < envs.size(); i++) {
    masses[i] = envs[i]->getRealEnvPtr()->getMonoNeutralMass() + mass_constant::getProtonMass();
  }
  return masses;
}

std::vector<double> getIntensitySums(const MatchEnvPtrVec &envs) {
  std::vector<double> intensity_sums(envs.size());
  for (size_t i = 0; i < envs.size(); i++) {
    intensity_sums[i] = envs[i]->getTheoEnvPtr()->compIntensitySum();
  }
  return intensity_sums;
}

void assignIntensity(PeakPtrVec &ms, MatchEnvPtrVec &envs) {
  size_t peak_num = ms.size();
  std::vector<double> intensity_sums_(peak_num, 0);
  for (size_t i = 0; i < envs.size(); i++) {
    MatchEnvPtr env = envs[i];
    for (int j = 0; j < env->getRealEnvPtr()->getPeakNum(); j++) {
      int peak = env->getRealEnvPtr()->getPeakIdx(j);
      if (peak >= 0) {
        intensity_sums_[peak] += env->getTheoEnvPtr()->getIntensity(j);
      }
    }
  }
  for (size_t i = 0; i < envs.size(); i++) {
    MatchEnvPtr env = envs[i];
    for (int j = 0; j < env->getRealEnvPtr()->getPeakNum(); j++) {
      int peak = env->getRealEnvPtr()->getPeakIdx(j);
      if (peak >= 0) {
        EnvelopePtr real_env = env->getRealEnvPtr();
        double intensity = real_env->getIntensity(j) * env->getTheoEnvPtr()->getIntensity(j)
            / intensity_sums_[peak];
        real_env->setIntensity(j, intensity);
      }
    }
  }
}

PeakPtrVec  rmAnnoPeak(PeakPtrVec &ms, MatchEnvPtrVec &envs) {
  PeakPtrVec new_list = ms;
  int peak_num = new_list.size();
  std::vector<bool> is_keeps(peak_num, true);
  for (size_t i = 0; i < envs.size(); i++) {
    RealEnvPtr real_env_ptr = envs[i]->getRealEnvPtr();
    for (int j = 0; j < real_env_ptr->getPeakNum(); j++) {
      int peak_idx = real_env_ptr->getPeakIdx(j);
      if (peak_idx >= 0) {
        is_keeps[peak_idx] = false;
      }
    }
  }
  raw_ms_util::rmPeaks(new_list, is_keeps);
  return new_list;
}

MatchEnvPtrVec addLowMassPeak(MatchEnvPtrVec &envs, PeakPtrVec &ms, double tolerance) {
  std::vector<bool> is_uses(ms.size(), false);
  for (size_t i = 0; i < envs.size(); i++) {
    MatchEnvPtr env = envs[i];
    for (int j = 0; j < env->getRealEnvPtr()->getPeakNum(); j++) {
      int peak = env->getRealEnvPtr()->getPeakIdx(j);
      if (peak >= 0) {
        is_uses[peak] = true;
      }
    }
  }

  MatchEnvPtrVec low_mass_envs;
  for (size_t i = 0; i < is_uses.size(); i++) {
    if (!is_uses[i]) {
      low_mass_envs.push_back(getNewMatchEnv(ms, i, tolerance));
    }
  }

  MatchEnvPtrVec result;
  result.insert(std::end(result), std::begin(low_mass_envs), std::end(low_mass_envs));
  result.insert(std::end(result), std::begin(envs), std::end(envs));
  std::sort(result.begin(), result.end(), MatchEnv::cmpScoreDec);
  return result;
}

MatchEnvPtr getNewMatchEnv(PeakPtrVec &ms, int idx, double tolerance) {
  double mz = ms[idx]->getPosition();
  double inte = ms[idx]->getIntensity();
  EnvPeakPtrVec peaks;
  EnvPeakPtr peak_ptr = std::make_shared<EnvPeak>(mz, inte);
  peaks.push_back(peak_ptr);
  int ref_idx = 0;
  int charge = 1;
  EnvelopePtr theo_env = std::make_shared<Envelope>(ref_idx, charge, mz, peaks);
  RealEnvPtr real_env = std::make_shared<RealEnv>(ms, theo_env, tolerance, 0);
  int mass_group = 0;
  return std::make_shared<MatchEnv>(mass_group, theo_env, real_env);
}

MatchEnvPtrVec addMultipleMass(MatchEnvPtrVec &envs, MatchEnvPtr2D &candidates,
                               double multi_min_mass, int multi_min_charge, double min_ratio) {
  MatchEnvPtrVec mass_envs;
  for (size_t i = 0; i < envs.size(); i++) {
    // check if we can use another charge state
    int charge = envs[i]->getRealEnvPtr()->getCharge();
    int refer_peak = envs[i]->getRealEnvPtr()->getReferPeakIdx();
    MatchEnvPtrVec charge_envs(2, nullptr);
    // we use non-overlapping envelopes here
    charge_envs[0] = candidates[refer_peak][charge-1];
    double min_score = charge_envs[0]->getScore() * min_ratio;
    if (charge >= multi_min_charge) {
      double score_minus_one  = 0;
      if (charge > 1 && candidates[refer_peak][charge-2] != nullptr) {
        score_minus_one = candidates[refer_peak][charge-2]->getScore();
      }
      double score_plus_one = 0;
      if (charge < static_cast<int>(candidates[0].size()) && candidates[refer_peak][charge] != nullptr) {
        score_plus_one = candidates[refer_peak][charge]->getScore();
      }
      if (score_minus_one > score_plus_one && score_minus_one >= min_score) {
        charge_envs[1] = candidates[refer_peak][charge-2];
      } else if (score_plus_one >= min_score) {
        charge_envs[1] = candidates[refer_peak][charge];
      }
    }
    for (size_t j = 0; j < charge_envs.size(); j++) {
      if (charge_envs[j] == nullptr) {
        continue;
      }
      mass_envs.push_back(charge_envs[j]);
      double mono_mass = charge_envs[j]->getRealEnvPtr()->getMonoNeutralMass();
      RealEnvPtr real_env_ptr = charge_envs[j]->getRealEnvPtr();
      int refer_idx = real_env_ptr->getReferIdx();
      if (mono_mass >= multi_min_mass) {
        /* check left shift */
        if (refer_idx > 0) {
          int p = real_env_ptr->getPeakIdx(refer_idx-1);
          if (p >=0 && candidates[p][charge-1] != nullptr &&
              candidates[p][charge-1]->getScore() >= min_score) {
            mass_envs.push_back(candidates[p][charge-1]);
          }
        }
        if (refer_idx < real_env_ptr->getPeakNum() - 1) {
          int p = real_env_ptr->getPeakIdx(refer_idx+1);
          if (p >=0 && candidates[p][charge-1] != nullptr &&
              candidates[p][charge-1]->getScore() >= min_score) {
            mass_envs.push_back(candidates[p][charge-1]);
          }
        }
      }
    }
  }
  return mass_envs;
}

DeconvMsPtr getDeconvMsPtr(MsHeaderPtr header_ptr, MatchEnvPtrVec &envs) {
  DeconvPeakPtrVec peak_list;
  int sp_id = header_ptr->getId();
  for (size_t i = 0; i < envs.size(); i++) {
    EnvelopePtr theo_env = envs[i]->getTheoEnvPtr();
    RealEnvPtr real_env = envs[i]->getRealEnvPtr();
    double pos = real_env->getMonoNeutralMass();
    double inte = theo_env->compIntensitySum();
    int charge = theo_env->getCharge();
    double score = envs[i]->getScore();
    DeconvPeakPtr peak_ptr = std::make_shared<DeconvPeak>(sp_id, i, pos, inte, charge, score);
    peak_list.push_back(peak_ptr);
  }
  DeconvMsPtr ms_ptr = std::make_shared<DeconvMs>(header_ptr, peak_list);
  return ms_ptr;
}

}  // namespace match_env_util

}  // namespace toppic
