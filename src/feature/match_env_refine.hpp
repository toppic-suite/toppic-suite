#ifndef PROT_FEATURE_MATCH_ENV_REFINE_HPP_
#define PROT_FEATURE_MATCH_ENV_REFINE_HPP_

#include "spec/peak.hpp"
#include "feature/match_env.hpp"

namespace prot {

class MatchEnvRefine {
 public:
  static double max_distance_a_;
  static double max_distance_b_;
  static double best_ratio_;

  static void mzRefine(FeatureMngPtr mng_ptr, MatchEnvPtrVec &envs);

  static void mzRefine(FeatureMngPtr mng_ptr, MatchEnvPtr env);

  static double compEnvDist(EnvelopePtr real_env, EnvelopePtr theo_env);

  static double compDistWithNorm(std::vector<double> real, std::vector<double> theo);

  static std::vector<double> norm(std::vector<double> &obs, double ratio);

  static double compDist(std::vector<double> &norm, std::vector<double> &theo);
};

}

#endif
