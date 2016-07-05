#include <algorithm>

#include "feature/env_assign.hpp" 

namespace prot {

MatchEnvPtr2D EnvAssign::assignWinEnv(MatchEnvPtr2D &match_envs, DeconvDataPtr data_ptr,
                          int env_num_per_win) {
  int win_num = data_ptr->getWinNum();
  MatchEnvPtr2D env_list(win_num);
  // add matchenv to the list */
  for (size_t i = 0; i < match_envs.size(); i++) {
    // i is peak idx 
    int win_id = data_ptr->getWinId(i);
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
