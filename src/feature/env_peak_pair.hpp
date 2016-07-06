#ifndef PROT_FEATURE_ENV_PEAK_PAIR_HPP_
#define PROT_FEATURE_ENV_PEAK_PAIR_HPP_

#include <memory>
#include <vector>

#include "feature/match_env.hpp"

namespace prot {

class EnvPeakPair;
typedef std::shared_ptr<EnvPeakPair> EnvPeakPairPtr;

class EnvPeakPair {
 public:
  EnvPeakPair(MatchEnvPtr env_ptr, int pos_idx);

  EnvPeakPair(EnvPeakPairPtr pair_ptr);

double getTheoIntensity();

double getPeakScore(double intensity_sum, double tolerance);

  MatchEnvPtr getMatchEnvPtr() {return env_ptr_;}
  int getPosIdx() { return pos_idx_;}

 private:
  MatchEnvPtr env_ptr_;
  int pos_idx_;
};

typedef std::vector<EnvPeakPairPtr> EnvPeakPairPtrVec;
typedef std::vector<EnvPeakPairPtrVec> EnvPeakPairPtr2D;

}
#endif
