#include "feature/env_peak_pair.hpp"

namespace prot {

EnvPeakPair::EnvPeakPair(MatchEnvPtr env_ptr, int pos_idx) {
  env_ptr_ = env_ptr;
  pos_idx_ = pos_idx;
}

EnvPeakPair::EnvPeakPair(EnvPeakPairPtr pair_ptr) {
  env_ptr_ = pair_ptr->getMatchEnvPtr();
  pos_idx_ = pair_ptr->getPosIdx();
}

double EnvPeakPair::getTheoIntensity() {
  return env_ptr_->getTheoEnvPtr()->getIntensity(pos_idx_);

}

double EnvPeakPair::getPeakScore(double intensity_sum, double tolerance) {
  return env_ptr_->calcPeakScr(pos_idx_, intensity_sum, tolerance);
}

}
