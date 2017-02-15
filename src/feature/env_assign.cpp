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
