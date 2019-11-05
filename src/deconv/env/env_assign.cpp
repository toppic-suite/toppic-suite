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


#include <algorithm>

#include "deconv/env/env_assign.hpp" 

namespace toppic {

namespace env_assign {

MatchEnvPtr2D assignWinEnv(MatchEnvPtr2D &match_envs, int win_num, 
                           const std::vector<int> &win_ids, int env_num_per_win) {
  MatchEnvPtr2D env_list(win_num);
  // add matchenv to the list */
  for (size_t i = 0; i < match_envs.size(); i++) {
    // i is peak idx 
    int win_id = win_ids[i];
    for (size_t j = 0; j < match_envs[i].size(); j++) {
      if (match_envs[i][j] != nullptr) {
        env_list[win_id].push_back(match_envs[i][j]);
      }
    }
  }
  // sort the matched envelopes and keep the best 
  for (int i = 0; i < win_num; i++) {
    std::sort(env_list[i].begin(), env_list[i].end(), MatchEnv::cmpScoreDec); 
    if ((int)env_list[i].size() > env_num_per_win) {
      env_list[i].resize(env_num_per_win);
    }
  }
  return env_list;
}

}

}
