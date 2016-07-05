#ifndef PROT_FEATURE_ENV_ASSIGN_HPP_
#define PROT_FEATURE_ENV_ASSIGN_HPP_

#include "feature/match_env.hpp"
#include "feature/deconv_data.hpp"

namespace prot {

class EnvAssign {
 public:
  static MatchEnvPtr2D assignWinEnv(MatchEnvPtr2D &match_envs, DeconvDataPtr data_ptr,
                                    int env_num_per_win);
};

}

#endif
