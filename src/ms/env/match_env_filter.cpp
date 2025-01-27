//Copyright (c) 2014 - 2025, The Trustees of Indiana University.
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
#include "ms/env/match_env_filter.hpp" 

namespace toppic {

namespace match_env_filter {

MatchEnvPtrVec filterByEnvCnnScore(MatchEnvPtrVec &ori_envs, double cutoff) {
  MatchEnvPtrVec result;
  for (size_t i = 0; i < ori_envs.size(); i++) {
    if (ori_envs[i]->getEnvcnnScore() >= cutoff) {
      result.push_back(ori_envs[i]);
    }
  }
  return result;
}

MatchEnvPtrVec filterByAANum(MatchEnvPtrVec &ori_envs, double prec_mass,
                             bool use_msdeconv, EnvParaPtr env_para_ptr) {
  MatchEnvPtrVec low_mass_envs;
  MatchEnvPtrVec high_mass_envs;
  if (use_msdeconv) {
    std::sort(ori_envs.begin(), ori_envs.end(), MatchEnv::cmpMsdeconvScoreDec);
  }
  else {
    std::sort(ori_envs.begin(), ori_envs.end(), MatchEnv::cmpEnvcnnScoreDec);
  }
  int low_mass_num = env_para_ptr->compLowMassNum(); 
  int high_mass_num = env_para_ptr->compHighMassNum(prec_mass);
  for (size_t i = 0; i < ori_envs.size(); i++) {
    if (ori_envs[i]->getExpEnvPtr()->getMonoNeutralMass() <= env_para_ptr->low_high_dividor_) {
      if ((int)low_mass_envs.size() < low_mass_num) {
        low_mass_envs.push_back(ori_envs[i]);
      }
    } else {
      if ((int)high_mass_envs.size() < high_mass_num) {
        high_mass_envs.push_back(ori_envs[i]);
      }
    }
  }
  MatchEnvPtrVec result;
  result.insert(std::end(result), std::begin(low_mass_envs), std::end(low_mass_envs));
  result.insert(std::end(result), std::begin(high_mass_envs), std::end(high_mass_envs));
  if (use_msdeconv) {
    std::sort(result.begin(), result.end(), MatchEnv::cmpMsdeconvScoreDec);
  }
  else {
    std::sort(result.begin(), result.end(), MatchEnv::cmpEnvcnnScoreDec);
  }
  return result;
}
}

}
