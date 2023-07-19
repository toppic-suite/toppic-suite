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

#include <cmath>
#include <algorithm>
#include <fstream>

#include "common/util/logger.hpp"
#include "ms/env/match_env.hpp"
#include "ms/env/env_detect.hpp"
#include "ms/env/env_filter.hpp"
#include "topfd/deconv/deconv_prec_win.hpp"

namespace toppic {

namespace deconv_prec_win {

struct PeakIntv {
  int bgn;
  int end;
};

const int PRECURSOR_TOP_PEAK_NUM = 20;

EnvParaPtr initMngPtr() {
  EnvParaPtr env_para_ptr = std::make_shared<EnvPara>();
  env_para_ptr->max_miss_peak_num_ = 3;
  std::vector<int> num = {1, 1, 1};
  env_para_ptr->min_match_peak_num_ = num;
  env_para_ptr->min_consecutive_peak_num_ = num;
  return env_para_ptr;
}

PeakIntv initPeakIntv(EnvParaPtr env_para_ptr, PeakPtrVec &peak_list, 
                      double bgn_mz, double end_mz) {
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
  if (min_dist < 0.01) {min_dist = 0.01;}
  
  int max_charge = static_cast<int>(std::round(1.0 / min_dist));

  if (max_charge > argu_max_charge) {
    max_charge = argu_max_charge;
  }
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
                           double argu_max_mass, int max_charge, 
                           double min_inte, double min_ref_inte) {
  MatchEnvPtr2D result;
  for (int idx = peak_intv.bgn; idx <= peak_intv.end; idx++) {
    MatchEnvPtrVec env_ptrs(max_charge);
    if (peak_list[idx]->getIntensity() >= min_inte) {
      for (int charge = 1; charge <= max_charge; charge++) {
        double max_mass = peak_list[idx]->getPosition() * charge + 1;
        MatchEnvPtr env_ptr;
        if (max_mass > argu_max_mass) {
          max_mass = argu_max_mass;
        } else {
          env_ptr  = env_detect::detectEnvByRefPeak(peak_list, idx, charge, 
                                                    max_mass, min_inte, min_ref_inte, 
                                                    env_para_ptr);
        }
        if (env_ptr != nullptr) {
          if (!env_filter::checkRealEnvValid(env_ptr, env_para_ptr)) {
            env_ptr = nullptr;
          } else {
            env_ptr->compMsdeconvScr(env_para_ptr);
          }
        }
        env_ptrs[charge-1] = env_ptr;
      }
    }
    result.push_back(env_ptrs);
  }
  return result;
}

MatchEnvPtr findBest(MatchEnvPtr2D &env_ptrs, double max_mass) {
  MatchEnvPtr best_env = nullptr;
  double best_score = -1;
  for (size_t i = 0; i < env_ptrs.size(); i++) {
    for (size_t j = 0; j < env_ptrs[i].size(); j++)  {
      if (env_ptrs[i][j] != nullptr && 
          env_ptrs[i][j]->getTheoEnvPtr()->getMonoNeutralMass() <= max_mass &&
          env_ptrs[i][j]->getMsdeconvScore() > best_score) {
        best_score = env_ptrs[i][j]->getMsdeconvScore();
        best_env = env_ptrs[i][j];
      }
    }
  }
  return best_env;
}

MatchEnvPtr deconv(double prec_win_begin, double prec_win_end, PeakPtrVec &peak_list,
                   double argu_max_mass, int argu_max_charge) {
  if (prec_win_begin <= 0) {
    return nullptr;
  }
  EnvParaPtr env_para_ptr = initMngPtr();
  PeakIntv peak_intv = initPeakIntv(env_para_ptr, peak_list, prec_win_begin, prec_win_end);
  int peak_num = initPeakNum(peak_intv);
  if (peak_num  == 0) {
    return nullptr;
  }
  int max_charge = initMaxChrg(peak_list, peak_intv, argu_max_charge);
  double min_inte = 0;
  double min_ref_inte = 0;
  LOG_DEBUG("Calcate match envelopes...");
  MatchEnvPtr2D match_envs = initMatchEnv(env_para_ptr, peak_list, peak_intv,
                                          peak_num, argu_max_mass, max_charge, min_inte, 
                                          min_ref_inte);
  LOG_DEBUG("Do filtering...");
  MatchEnvPtr env_ptr = findBest(match_envs, argu_max_mass);
  return env_ptr;
}

MatchEnvPtr deconvPrecWinForOneMs(MzmlMsPtr ms_one, MzmlMsPtr ms_two, 
                                  double max_mass, int max_charge) {
  MsHeaderPtr ms_two_header = ms_two->getMsHeaderPtr();
  double prec_win_begin = ms_two_header->getPrecWinBegin();
  double prec_win_end = ms_two_header->getPrecWinEnd();

  PeakPtrVec peak_list = ms_one->getPeakPtrVec();
  MatchEnvPtr match_env_ptr = deconv_prec_win::deconv(prec_win_begin, prec_win_end, peak_list,  
                                                      max_mass, max_charge);
  if (match_env_ptr == nullptr) {
    LOG_INFO("No precursor isotopic envelope is found for scan " << ms_two_header->getFirstScanNum());
  }
  return match_env_ptr;
}

MatchEnvPtrVec deconvPrecWinForMsGroup(MzmlMsGroupPtr ms_group_ptr, 
                                       double max_mass, int max_charge) {
  MzmlMsPtr ms_one_ptr = ms_group_ptr->getMsOnePtr();
  MzmlMsPtrVec ms_two_ptr_vec = ms_group_ptr->getMsTwoPtrVec();
  MatchEnvPtrVec result_envs;
  for (size_t i = 0; i < ms_two_ptr_vec.size(); i++) {
    MzmlMsPtr ms_two_ptr = ms_two_ptr_vec[i];
    MatchEnvPtr match_env_ptr = deconvPrecWinForOneMs(ms_one_ptr, ms_two_ptr, 
                                                      max_mass, max_charge);
    if (match_env_ptr != nullptr) {
      result_envs.push_back(match_env_ptr);
    }
  }
  return result_envs;
}

}

}  // namespace toppic
