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


#ifndef TOPPIC_TOPFD_DP_CO_TABLE_HPP_
#define TOPPIC_TOPFD_DP_CO_TABLE_HPP_

#include <vector>

#include "ms/env/match_env.hpp"

namespace toppic {

class CoTable {
 public:
  static int cntTableSize(MatchEnvPtr2D &win_envs, int win, int pnt);

  static bool checkCoexist(MatchEnvPtr env_a, MatchEnvPtr env_b, int pos_a,
                                    int pos_b, double tolerance);

  static bool checkCoexist(MatchEnvPtr env_a, MatchEnvPtr env_b,
                                    double tolerance);

  static void compTableEntry(MatchEnvPtrVec &env_list, MatchEnvPtr2D &win_envs,
                             std::vector<bool> &rows, int id, int win,
                             double tolerance);

  static std::vector<std::vector<bool>> initCoexistTable(MatchEnvPtr2D &win_envs, double tolerance);
};

}  // namespace toppic
#endif
