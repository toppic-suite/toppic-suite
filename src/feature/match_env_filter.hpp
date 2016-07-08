#ifndef PROT_FEATURE_MATCH_ENV_FILTER_HPP_
#define PROT_FEATURE_MATCH_ENV_FILTER_HPP_

#include "spec/peak.hpp"
#include "feature/match_env.hpp"

namespace prot {

class MatchEnvFilter {
 public:
  static MatchEnvPtrVec filter(MatchEnvPtrVec &ori_envs, double prec_mass, FeatureMngPtr mng_ptr);
};

}

#endif
