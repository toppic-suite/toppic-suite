// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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
#include <cstddef>
#include <limits>

#include "base/logger.hpp"
#include "feature/raw_ms_util.hpp" 
#include "feature/match_env_util.hpp" 

namespace prot {

std::vector<double> MatchEnvUtil::getMassList(const MatchEnvPtrVec &envs) {
  std::vector<double> masses(envs.size());
  for (size_t i = 0; i < envs.size(); i++) {
    masses[i] = envs[i]->getRealEnvPtr()->getMonoMass();
  }
  return masses;
}

std::vector<int> MatchEnvUtil::getChargeList(const MatchEnvPtrVec &envs) {
  std::vector<int> charges(envs.size());
  for (size_t i = 0; i < envs.size(); i++) {
    charges[i] = envs[i]->getRealEnvPtr()->getCharge();
  }
  return charges;
}

std::vector<double> MatchEnvUtil::getChargeOneMassList(const MatchEnvPtrVec &envs) {
  std::vector<double> masses(envs.size());
  for (size_t i = 0; i < envs.size(); i++) {
    masses[i] = envs[i]->getRealEnvPtr()->getMonoMass() + MassConstant::getProtonMass();
  }
  return masses;
}

std::vector<double> MatchEnvUtil::getIntensitySums(const MatchEnvPtrVec &envs) {
  std::vector<double> intensity_sums(envs.size());
  for (size_t i = 0; i < envs.size(); i++) {
    intensity_sums[i] = envs[i]->getTheoEnvPtr()->compIntensitySum();
  }
  return intensity_sums;
}

void MatchEnvUtil::assignIntensity(PeakPtrVec &ms, MatchEnvPtrVec &envs) {
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

PeakPtrVec  MatchEnvUtil::rmAnnoPeak(PeakPtrVec &ms, MatchEnvPtrVec &envs) {
  PeakPtrVec new_list = ms;
  int peak_num = new_list.size();
  std::vector<bool> is_keeps(peak_num, true);
  for (size_t i = 0; i < envs.size(); i++) {
    std::vector<int> peaks = envs[i]->getRealEnvPtr()->getPeakIdxList();
    for (size_t j = 0; j < peaks.size(); j++) {
      if (peaks[j] >= 0) {
        is_keeps[peaks[j]] = false;
      }
    }
  }
  RawMsUtil::rmPeaks(new_list, is_keeps);
  return new_list;
}


MatchEnvPtrVec MatchEnvUtil::addLowMassPeak(MatchEnvPtrVec &envs, PeakPtrVec &ms, 
                                            double tolerance) {
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

MatchEnvPtr MatchEnvUtil::getNewMatchEnv(PeakPtrVec &ms, int idx, double tolerance) {
  std::vector<double> mzs; 
  std::vector<double> intensities;
  mzs.push_back(ms[idx]->getPosition());
  intensities.push_back(ms[idx]->getIntensity());
  LOG_DEBUG(mzs[0] << "\t" << intensities[0]);
  EnvelopePtr theo_env = std::make_shared<Envelope>(0, 1, mzs[0], mzs, intensities);
  RealEnvPtr real_env = std::make_shared<RealEnv>(ms, theo_env, tolerance, 0);
  int mass_group = 0;
  return std::make_shared<MatchEnv>(mass_group, theo_env, real_env);
}

MatchEnvPtrVec MatchEnvUtil::addMultipleMass(MatchEnvPtrVec &envs, MatchEnvPtr2D &candidates,
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
      if (charge < (int)candidates[0].size() && candidates[refer_peak][charge] != nullptr) {
        score_plus_one = candidates[refer_peak][charge]->getScore();
      }
      if (score_minus_one > score_plus_one && score_minus_one >= min_score) {
        charge_envs[1] = candidates[refer_peak][charge-2];
      }
      else if (score_plus_one >= min_score) {
        charge_envs[1] = candidates[refer_peak][charge];
      }
    }
    for (size_t j = 0; j < charge_envs.size(); j++) {
      if (charge_envs[j] == nullptr) {
        continue;
      }
      mass_envs.push_back(charge_envs[j]);
      double mono_mass = charge_envs[j]->getRealEnvPtr()->getMonoMass();
      std::vector<int> peaks = charge_envs[j]->getRealEnvPtr()->getPeakIdxList();
      int refer_idx = charge_envs[j]->getRealEnvPtr()->getReferIdx();
      if (mono_mass >= multi_min_mass) {
        /* check left shift */
        if (refer_idx > 0) {
          int p = peaks[refer_idx-1];
          if (p >=0 && candidates[p][charge-1] != nullptr && 
              candidates[p][charge-1]->getScore() >= min_score) {
            mass_envs.push_back(candidates[p][charge-1]);
          }
        }
        if (refer_idx < (int)peaks.size() - 1) {
          int p = peaks[refer_idx+1];
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

}  // namespace prot
