#ifndef PROT_FEATURE_ENV_FILTER_HPP_
#define PROT_FEATURE_ENV_FILTER_HPP_

#include "feature/match_env.hpp"
#include "feature/deconv_data.hpp"

namespace prot {

class EnvFilter {
 public:
  static void filter(MatchEnvPtr2D &match_envs, DeconvDataPtr data_ptr,
                     FeatureMngPtr mng_ptr);
  
  static void multipleMassFilter(MatchEnvPtr2D &match_env, DeconvDataPtr data_ptr,
                                 FeatureMngPtr mng_ptr);

  static bool testRealEnvValid(MatchEnvPtr env, FeatureMngPtr mng_ptr);
};

}

#endif
