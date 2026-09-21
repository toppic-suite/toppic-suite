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

#include "ms/env/env_detect.hpp"

#include <climits>
#include <cmath>
#include <memory>
#include <vector>

#include "common/base/mass_constant.hpp"
#include "common/util/logger.hpp"
#include "ms/env/env_base.hpp"
#include "ms/spec/peak_list_util.hpp"

namespace toppic {

namespace env_detect {

namespace {

// compute the intensity ratio based on the top three peaks
double calcInteRatio(const PeakPtrVec& peak_list, EnvPtr theo_env_ptr,
                     double tolerance) {
  double theo_sum = 0;
  double obs_sum = 0;
  int refer_idx = theo_env_ptr->getReferIdx();
  double mz = theo_env_ptr->getMz(refer_idx);
  int peak_idx = peak_list_util::getNearPeakIdx(peak_list, mz, tolerance);
  if (peak_idx >= 0) {
    theo_sum += theo_env_ptr->getInte(refer_idx);
    obs_sum += peak_list[peak_idx]->getIntensity();
  }
  if (refer_idx - 1 >= 0) {
    theo_sum += theo_env_ptr->getInte(refer_idx - 1);
    mz = theo_env_ptr->getMz(refer_idx - 1);
    peak_idx = peak_list_util::getNearPeakIdx(peak_list, mz, tolerance);
    if (peak_idx >= 0) {
      obs_sum += peak_list[peak_idx]->getIntensity();
    }
  }
  if (refer_idx + 1 < theo_env_ptr->getPeakNum()) {
    theo_sum += theo_env_ptr->getInte(refer_idx + 1);
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

MatchEnvPtr compEnv(const PeakPtrVec& peak_list, EnvPtr theo_env_ptr,
                    const EnvParaPtr& env_para_ptr, int mass_group,
                    double min_inte, double bound_inte) {
  // scale theoretical distribution
  int charge = theo_env_ptr->getCharge();
  double ratio = calcInteRatio(peak_list, theo_env_ptr,
                               env_para_ptr->getMzTolerance(charge));
  theo_env_ptr->changeIntensity(ratio);

  // get the highest 85%--95% peaks above the noise level
  double percentage = env_para_ptr->getPercentBound(mass_group);
  theo_env_ptr = theo_env_ptr->getSubEnv(percentage, bound_inte,
                                         env_para_ptr->max_back_peak_num_,
                                         env_para_ptr->max_forw_peak_num_);
  // get real envelope
  ExpEnvPtr real_env_ptr = std::make_shared<ExpEnv>(
      peak_list, theo_env_ptr, env_para_ptr->getMzTolerance(charge), min_inte);
  MatchEnvPtr match_env_ptr =
      std::make_shared<MatchEnv>(mass_group, theo_env_ptr, real_env_ptr);
  return match_env_ptr;
}

MatchEnvPtr compFullEnv(const PeakPtrVec& peak_list, EnvPtr theo_env_ptr,
                        const EnvParaPtr& env_para_ptr, int mass_group,
                        double min_inte, double bound_inte) {
  // scale theoretical distribution
  int charge = theo_env_ptr->getCharge();
  double ratio = calcInteRatio(peak_list, theo_env_ptr,
                               env_para_ptr->getMzTolerance(charge));
  theo_env_ptr->changeIntensity(ratio);

  double percentage = 1.0;
  double max_back = INT_MAX;
  double max_forw = INT_MAX;
  theo_env_ptr =
      theo_env_ptr->getSubEnv(percentage, bound_inte, max_back, max_forw);
  // get real envelope
  ExpEnvPtr real_env_ptr = std::make_shared<ExpEnv>(
      peak_list, theo_env_ptr, env_para_ptr->getMzTolerance(charge), min_inte);
  MatchEnvPtr match_env_ptr =
      std::make_shared<MatchEnv>(mass_group, theo_env_ptr, real_env_ptr);
  return match_env_ptr;
}

// Exact pre-check of env_filter::checkRealEnvValid for a candidate of
// getCandidateEnv, computed without building the theoretical envelope, the
// experimental envelope and the MatchEnv: the same shifted positions
// (Env::convertToTheo), the same intensity ratio (calcInteRatio), the same
// bound of peaks above the noise (Env::calcBound with percentage 1 and
// unlimited peak numbers, as in compFullEnv), the same nearest-peak mapping
// and duplicate-match rule (ExpEnv::mapPeakList / remvDuplMatch) and the
// same matched, missing and consecutive counts. Returns false only when the
// filter would remove the envelope, so the candidate lists that reach the
// later filters and the DP are unchanged. The filter keeps about 1% of the
// candidates (toppic_test experiment 36), so this skips the allocations of
// three envelopes with their peaks for the other 99%.
bool passesRealEnvCheck(const PeakPtrVec& peak_list, const EnvPtr& ref_env_ptr,
                        double refer_mz, int charge, int mass_group,
                        double min_inte, double bound_inte,
                        const EnvParaPtr& env_para_ptr) {
  int ori_charge = 1;
  int n = ref_env_ptr->getPeakNum();
  int refer_idx = ref_env_ptr->getReferIdx();
  // distrToTheoRef / convertToTheo
  double mass_diff =
      refer_mz * charge - ref_env_ptr->getReferMz() * ori_charge;
  std::vector<double> mz(n);
  std::vector<double> inte(n);
  for (int i = 0; i < n; i++) {
    mz[i] = (ref_env_ptr->getMz(i) + mass_diff) / charge;
    inte[i] = ref_env_ptr->getInte(i);
  }
  double tolerance = env_para_ptr->getMzTolerance(charge);
  // calcInteRatio
  double theo_sum = 0;
  double obs_sum = 0;
  int peak_idx = peak_list_util::getNearPeakIdx(peak_list, mz[refer_idx],
                                                tolerance);
  if (peak_idx >= 0) {
    theo_sum += inte[refer_idx];
    obs_sum += peak_list[peak_idx]->getIntensity();
  }
  if (refer_idx - 1 >= 0) {
    theo_sum += inte[refer_idx - 1];
    peak_idx = peak_list_util::getNearPeakIdx(peak_list, mz[refer_idx - 1],
                                              tolerance);
    if (peak_idx >= 0) {
      obs_sum += peak_list[peak_idx]->getIntensity();
    }
  }
  if (refer_idx + 1 < n) {
    theo_sum += inte[refer_idx + 1];
    peak_idx = peak_list_util::getNearPeakIdx(peak_list, mz[refer_idx + 1],
                                              tolerance);
    if (peak_idx >= 0) {
      obs_sum += peak_list[peak_idx]->getIntensity();
    }
  }
  double ratio = (theo_sum == 0) ? 1.0 : obs_sum / theo_sum;
  // changeIntensity
  for (int i = 0; i < n; i++) {
    inte[i] = inte[i] * ratio;
  }
  // calcBound(1.0, bound_inte, INT_MAX, INT_MAX)
  double percent_bound = 1.0;
  double sum = inte[refer_idx];
  int forw_idx = refer_idx + 1;
  int back_idx = refer_idx - 1;
  double intensity_sum = 0;
  for (int i = 0; i < n; i++) {
    intensity_sum = intensity_sum + inte[i];
  }
  while (sum / intensity_sum < percent_bound) {
    if (back_idx == -1 && forw_idx == n) {
      break;
    }
    double forw_inte = (forw_idx == n) ? 0 : inte[forw_idx];
    double back_inte = (back_idx == -1) ? 0 : inte[back_idx];
    if (back_inte < bound_inte && forw_inte < bound_inte) {
      break;
    }
    if (back_inte > forw_inte) {
      sum += back_inte;
      back_idx--;
    } else {
      sum += forw_inte;
      forw_idx++;
    }
  }
  int n_back = refer_idx - back_idx - 1;
  int n_forw = forw_idx - refer_idx - 1;
  int bgn = refer_idx - n_back;
  int m = n_back + n_forw + 1;
  // ExpEnv::mapPeakList
  std::vector<int> idx(m);
  std::vector<double> pos(m);
  for (int i = 0; i < m; i++) {
    double theo_mz = mz[bgn + i];
    int p = peak_list_util::getNearPeakIdx(peak_list, theo_mz, tolerance);
    if (p >= 0 && peak_list[p]->getIntensity() >= min_inte) {
      idx[i] = p;
      pos[i] = peak_list[p]->getPosition();
    } else {
      idx[i] = EnvPeak::getNonExistPeakIdx();
      pos[i] = theo_mz;
    }
  }
  // ExpEnv::remvDuplMatch
  int non_exist = EnvPeak::getNonExistPeakIdx();
  for (int i = 0; i < m - 1; i++) {
    if (idx[i] != non_exist && idx[i] == idx[i + 1]) {
      if (std::abs(mz[bgn + i] - pos[i]) <
          std::abs(mz[bgn + i + 1] - pos[i + 1])) {
        idx[i + 1] = non_exist;
      } else {
        idx[i] = non_exist;
      }
    }
  }
  // cntMissPeakNum, cntMaxConsPeakNum
  int miss_num = 0;
  int max_cons = 0;
  int cons = 0;
  for (int i = 0; i < m; i++) {
    if (idx[i] == non_exist) {
      miss_num++;
      cons = 0;
    } else {
      cons++;
      if (cons > max_cons) {
        max_cons = cons;
      }
    }
  }
  // checkRealEnvValid
  if (m - miss_num < env_para_ptr->min_match_peak_num_[mass_group]) {
    return false;
  }
  if (miss_num > env_para_ptr->max_miss_peak_num_) {
    return false;
  }
  if (env_para_ptr->check_consecutive_peak_num_) {
    if (max_cons < env_para_ptr->compMinConsPeakNum(m, mass_group)) {
      return false;
    }
  }
  return true;
}

}  // namespace

// detect a MatchEnv
MatchEnvPtr detectEnvByRefPeak(const PeakPtrVec& peak_list, int ref_peak,
                               int charge, double max_mass, double min_inte,
                               double min_ref_inte, double bound_inte,
                               const EnvParaPtr& env_para_ptr, bool is_full,
                               bool check_real_env) {
  double refer_mz = peak_list[ref_peak]->getPosition();
  // check if the mass is greater than the precursor mass
  double ref_mass = peak_util::compPeakNeutralMass(refer_mz, charge);
  if (ref_mass >= max_mass || ref_mass < env_para_ptr->min_mass_) {
    return nullptr;
  }

  // get a reference distribution based on the reference (highest intensity)
  // mass
  EnvPtr ref_env_ptr = EnvBase::getEnvByRefMass(ref_mass);
  if (ref_env_ptr == nullptr) {
    LOG_WARN("reference envelope is null");
    return nullptr;
  }

  if (peak_list[ref_peak]->getIntensity() < min_ref_inte) {
    return nullptr;
  }
  // the m/z of the reference peak of the theoretical distribution shifted
  // to refer_mz and charge (Env::convertToTheo applied to the reference
  // peak: (ref_mz + (refer_mz * charge - ref_mz)) / charge)
  int ori_charge = 1;
  double mass_diff = refer_mz * charge - ref_env_ptr->getReferMz() * ori_charge;
  double theo_refer_mz = (ref_env_ptr->getReferMz() + mass_diff) / charge;
  int peak_idx = peak_list_util::getNearPeakIdx(
      peak_list, theo_refer_mz, env_para_ptr->getMzTolerance(charge));
  if (peak_idx < 0 || peak_list[peak_idx]->getIntensity() < min_ref_inte) {
    return nullptr;
  }
  int mass_group = env_para_ptr->getMassGroup(ref_mass);
  if (is_full && check_real_env &&
      !passesRealEnvCheck(peak_list, ref_env_ptr, refer_mz, charge,
                          mass_group, min_inte, bound_inte, env_para_ptr)) {
    return nullptr;
  }
  // convert the distribution to a theoretical distribution
  // based on the refer mz and charge state
  EnvPtr theo_env_ptr = ref_env_ptr->distrToTheoRef(refer_mz, charge);
  if (is_full) {
    return compFullEnv(peak_list, theo_env_ptr, env_para_ptr, mass_group,
                       min_inte, bound_inte);
  } else {
    return compEnv(peak_list, theo_env_ptr, env_para_ptr, mass_group, min_inte,
                   bound_inte);
  }
}

MatchEnvPtr detectEnvByMonoMass(const PeakPtrVec& peak_list, double mono_mass,
                                int charge, double min_inte, double bound_inte,
                                const EnvParaPtr& env_para_ptr, bool is_full) {
  if (mono_mass < env_para_ptr->min_mass_) {
    return nullptr;
  }

  // get a reference distribution based on the base mass
  EnvPtr ref_env_ptr = EnvBase::getEnvByMonoMass(mono_mass);
  if (ref_env_ptr == nullptr) {
    LOG_WARN("Reference envelope is null!");
    return nullptr;
  }

  // convert the reference distribution to a theoretical distribution
  // based on the mono mz and charge state
  double mono_mz = mono_mass / charge + mass_constant::getProtonMass();
  EnvPtr theo_env_ptr = ref_env_ptr->distrToTheoMono(mono_mz, charge);

  int mass_group = env_para_ptr->getMassGroup(mono_mass);
  if (is_full) {
    return compFullEnv(peak_list, theo_env_ptr, env_para_ptr, mass_group,
                       min_inte, bound_inte);
  } else {
    return compEnv(peak_list, theo_env_ptr, env_para_ptr, mass_group, min_inte,
                   bound_inte);
  }
}

MatchEnvPtr2D getCandidateEnv(const PeakPtrVec& peak_list, int max_charge,
                              double max_mass, double min_inte,
                              double min_ref_inte, double bound_inte,
                              const EnvParaPtr& env_para_ptr,
                              bool check_real_env) {
  bool is_full = true;
  int peak_num = peak_list.size();
  MatchEnvPtr2D match_envs(peak_num);
  for (int idx = 0; idx < peak_num; idx++) {
    match_envs[idx].resize(max_charge);
    for (int charge = 1; charge <= max_charge; charge++) {
      // with check_real_env, the candidates that env_filter::filter's
      // first test would remove are skipped before they are built
      match_envs[idx][charge - 1] = detectEnvByRefPeak(
          peak_list, idx, charge, max_mass, min_inte, min_ref_inte,
          bound_inte, env_para_ptr, is_full, check_real_env);
    }
  }
  return match_envs;
}

}  // namespace env_detect

}  // namespace toppic
