//Copyright (c) 2014 - 2023, The Trustees of Indiana University.
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
#include <algorithm>

#include "common/util/logger.hpp"
#include "common/base/mass_constant.hpp"

#include "topfd/ecscore/env/seed_env_util.hpp"
#include "topfd/ecscore/env/ms_map_env_util.hpp"
#include "topfd/ecscore/env_set/env_set.hpp"
#include "topfd/ecscore/env_set/env_set_util.hpp"

namespace toppic {

namespace env_set_util {

void removeNonMatchEnvs(MsMapEnvPtrVec &env_list, int refer_idx,
                        int min_match_peak_num) {
  int idx = env_list.size() - 1;
  while (idx >= 0) {
    MsMapEnvPtr env = env_list[idx];
    if (env->getTopThreeMatchNum(refer_idx) < min_match_peak_num)
      env_list.erase(env_list.begin() + idx);
    else
      return;
    idx = idx - 1;
  }
}

EnvSetPtr searchEnvSet(MsMapPtr ms_map_ptr, SeedEnvPtr seed_ptr,
                      EcscoreParaPtr para_ptr, double sn_ratio) { 
  int start_spec_id = 0;
  int end_spec_id = ms_map_ptr->getRowNum() - 1;
  return searchEnvSet(ms_map_ptr, seed_ptr, 
                      start_spec_id, end_spec_id,
                      para_ptr, sn_ratio);
}

EnvSetPtr searchEnvSet(MsMapPtr ms_map_ptr, SeedEnvPtr seed_ptr,
                       int start_spec_id, int end_spec_id,
                       EcscoreParaPtr para_ptr, double sn_ratio) {
  double mass_tole = para_ptr->getMassTole();
  int refer_peak_idx = seed_ptr->getReferIdx();
  double min_inte = ms_map_ptr->getBaseInte() * sn_ratio;
   // search backward
  MsMapEnvPtrVec back_env_list;
  int miss_num = 0;
  int spec_id = seed_ptr->getSpecId();
  while (spec_id >= start_spec_id) {
    MsMapEnvPtr ms_map_env_ptr = ms_map_env_util::getMatchMsMapEnv(ms_map_ptr, seed_ptr,
                                                                   spec_id, mass_tole, min_inte);
    back_env_list.push_back(ms_map_env_ptr);
    if (ms_map_env_ptr->getTopThreeMatchNum(refer_peak_idx)
       < para_ptr->min_match_peak_) {
      miss_num = miss_num + 1;
    }
    else {
      miss_num = 0;
    }
    if (miss_num >= para_ptr->max_miss_env_) {
      break;
    }
    spec_id = spec_id - 1;
  }
  removeNonMatchEnvs(back_env_list, refer_peak_idx, para_ptr->min_match_peak_);

  // search forward
  MsMapEnvPtrVec forw_env_list;
  spec_id = seed_ptr->getSpecId() + 1;
  miss_num = 0;
  while (spec_id <= end_spec_id) {
    MsMapEnvPtr ms_map_env_ptr = ms_map_env_util::getMatchMsMapEnv(ms_map_ptr, seed_ptr,
                                                                   spec_id, mass_tole, min_inte);
    forw_env_list.push_back(ms_map_env_ptr);
    if (ms_map_env_ptr->getTopThreeMatchNum(refer_peak_idx) < para_ptr->min_match_peak_) {
      miss_num = miss_num + 1;
    }
    else {
      miss_num = 0;
    }
    if (miss_num >= para_ptr->max_miss_env_) {
      break;
    }
    spec_id = spec_id + 1;
  }
  removeNonMatchEnvs(forw_env_list, refer_peak_idx, para_ptr->min_match_peak_);
  // merge results
  std::reverse(back_env_list.begin(), back_env_list.end());
  back_env_list.insert(back_env_list.end(), forw_env_list.begin(), forw_env_list.end());
  if (back_env_list.empty()) {
    return nullptr;
  }
  start_spec_id = back_env_list[0]->getSpecId();
  end_spec_id = back_env_list[back_env_list.size() - 1]->getSpecId();
  if ((end_spec_id - start_spec_id + 1) < para_ptr->min_match_env_) {
    return nullptr;
  }
  EnvSetPtr env_set_ptr = std::make_shared<EnvSet>(seed_ptr, back_env_list, 
                                                   start_spec_id, end_spec_id, 
                                                   min_inte);
  return env_set_ptr;
}

}

}
