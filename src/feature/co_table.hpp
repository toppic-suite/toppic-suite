#ifndef PROT_FEATURE_CO_TABLE_HPP_
#define PROT_FEATURE_CO_TABLE_HPP_

#include "feature/match_env.hpp"

namespace prot {

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

}
#endif
