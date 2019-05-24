//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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
#include <algorithm>

#include "common/util/logger.hpp"
#include "deconv/env/match_env.hpp"
#include "deconv/env/env_detect.hpp"
#include "deconv/env/env_filter.hpp"
#include "deconv/env/prec_env.hpp"

namespace toppic {

namespace prec_env {

struct PeakIntv {
  int bgn;
  int end;
};

const int PRECURSOR_TOP_PEAK_NUM = 20;

EnvParaPtr initMngPtr(double prec_win_size) {
  EnvParaPtr env_para_ptr = std::make_shared<EnvPara>();
  env_para_ptr->min_refer_inte_ = 0;
  env_para_ptr->min_inte_ = 0;
  env_para_ptr->max_miss_peak_num_ = 3;
  std::vector<int> num = {1, 1, 1};
  env_para_ptr->min_match_peak_num_ = num;
  env_para_ptr->min_consecutive_peak_num_ = num;
  env_para_ptr->prec_deconv_interval_ = prec_win_size;
  return env_para_ptr;
}

PeakIntv initPeakIntv(EnvParaPtr env_para_ptr, PeakPtrVec &peak_list, double prec_mz) {
  double base_mz  = prec_mz;
  double bgn_mz = base_mz - env_para_ptr->prec_deconv_interval_ / 2;
  double end_mz = base_mz + env_para_ptr->prec_deconv_interval_ / 2;
  PeakIntv peak_intv;
  peak_intv.bgn = peak_list.size();
  peak_intv.end = -1;
  for (size_t i = 0; i < peak_list.size(); i++) {
    if (peak_list[i]->getPosition() >= bgn_mz
        && peak_list[i]->getPosition() <= end_mz) {
      if (static_cast<int>(i) < peak_intv.bgn) {
        peak_intv.bgn = i;
      }
      if (static_cast<int>(i) > peak_intv.end) {
        peak_intv.end = i;
      }
    }
  }
  return peak_intv;
}

int initPeakNum(PeakIntv peak_intv) {
  int peak_num;
  if (peak_intv.end < peak_intv.bgn) {
    peak_num = 0;
  } else {
    peak_num = peak_intv.end - peak_intv.bgn + 1;
  }
  return peak_num;
}

int initMaxChrg(PeakPtrVec &peak_list, PeakIntv peak_intv, int argu_max_charge) {
  double min_dist = 1;
  for (int i = peak_intv.bgn - 1; i <= peak_intv.end; i++) {
    if (i < 0) {
      continue;
    }
    double cur_mz = peak_list[i]->getPosition();
    if (i + 1 >= static_cast<int>(peak_list.size())) {
      continue;
    }
    double next_mz = peak_list[i+1]->getPosition();
    double dist = next_mz - cur_mz;
    if (dist < min_dist) {
      min_dist = dist;
    }
  }
  int max_charge = static_cast<int>(std::round(1.0 / min_dist));
  if (max_charge > argu_max_charge) {
    max_charge = argu_max_charge;
  }
  // LOG_DEBUG("maximum charge: " << max_charge);
  return max_charge;
}

double initMinInte(PeakPtrVec &peak_list, 
                   PeakIntv peak_intv) { 
  int peak_num = peak_intv.end - peak_intv.bgn + 1;
  if (peak_num < PRECURSOR_TOP_PEAK_NUM) {
    return 0;
  }
  else {
    std::vector<PeakPtr>::const_iterator first = peak_list.begin() + peak_intv.bgn;
    std::vector<PeakPtr>::const_iterator last = peak_list.begin() + peak_intv.end + 1;
    PeakPtrVec new_peak_list(first, last);
    std::sort(new_peak_list.begin(), new_peak_list.end(), Peak::cmpInteDec);
    return new_peak_list[PRECURSOR_TOP_PEAK_NUM-1]->getIntensity();

  }
}

MatchEnvPtr2D initMatchEnv(EnvParaPtr env_para_ptr, PeakPtrVec &peak_list,
                           PeakIntv peak_intv, int peak_num, 
                           int max_charge, double min_inte) {
  MatchEnvPtr2D result;
  // LOG_DEBUG("peak intv bgn " << peak_intv.bgn << " end " << peak_intv.end);
  for (int idx = peak_intv.bgn; idx <= peak_intv.end; idx++) {
    MatchEnvPtrVec env_ptrs(max_charge);
    // LOG_DEBUG("idx " << idx << " mz " << peak_list[idx]->getPosition());
    if (peak_list[idx]->getIntensity() >= min_inte) {
      for (int charge = 1; charge <= max_charge; charge++) {
        double max_mass = peak_list[idx]->getPosition() * charge + 1;
        // LOG_DEBUG("max_mass " << max_mass);
        MatchEnvPtr env_ptr;
        if (max_mass > env_para_ptr->max_mass_) {
          max_mass = env_para_ptr->max_mass_;
        } else {
          env_ptr  = EnvDetect::detectEnv(peak_list, idx, charge, max_mass, env_para_ptr);
        }
        // LOG_DEBUG("env detection complete");
        if (env_ptr != nullptr) {
          // LOG_DEBUG("FILTER REAL ENVELOPE!!!");
          if (!EnvFilter::testRealEnvValid(env_ptr, env_para_ptr)) {
            env_ptr = nullptr;
          } else {
            env_ptr->compScr(env_para_ptr);
            // LOG_DEBUG("GOOD ENVELOPE FOUNDER!!!");
          }
        }
        env_ptrs[charge-1] = env_ptr;
      }
    }
    result.push_back(env_ptrs);
  }
  return result;
}

MatchEnvPtr findBest(MatchEnvPtr2D &env_ptrs) {
  MatchEnvPtr best_env = nullptr;
  double best_score = -1;
  for (size_t i = 0; i < env_ptrs.size(); i++) {
    for (size_t j = 0; j < env_ptrs[i].size(); j++)  {
      if (env_ptrs[i][j] != nullptr && env_ptrs[i][j]->getScore() > best_score) {
        best_score = env_ptrs[i][j]->getScore();
        best_env = env_ptrs[i][j];
      }
    }
  }
  return best_env;
}

MatchEnvPtr deconv(double prec_win_size, PeakPtrVec &peak_list,
                  double prec_mz, int prec_charge, int argu_max_charge) {
  LOG_DEBUG("Prec: " << prec_mz << " charge: " << prec_charge);
  if (prec_mz <= 0) {
    return nullptr;
  }
  EnvParaPtr env_para_ptr = initMngPtr(prec_win_size);
  PeakIntv peak_intv = initPeakIntv(env_para_ptr, peak_list, prec_mz);
  int peak_num = initPeakNum(peak_intv);
  if (peak_num  == 0) {
    return nullptr;
  }
  int max_charge = initMaxChrg(peak_list, peak_intv, argu_max_charge);
  //double min_inte = initMinInte(peak_list, peak_intv);
  double min_inte = 0;
  LOG_DEBUG("Calcate match envelopes...");
  MatchEnvPtr2D match_envs = initMatchEnv(env_para_ptr, peak_list, peak_intv,
                                          peak_num, max_charge, min_inte);
  LOG_DEBUG("Do filtering...");
  MatchEnvPtr env_ptr = findBest(match_envs);
  return env_ptr;
  /*
  if (env_ptr != nullptr) {
    return env_ptr;->getRealEnvPtr();
  } else {
    return nullptr;
  }
  */
}

}

}  // namespace toppic
