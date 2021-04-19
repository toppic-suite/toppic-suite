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

#ifndef TOPPIC_TOPFD_ENV_ENV_PEAK_PAIR_HPP_
#define TOPPIC_TOPFD_ENV_ENV_PEAK_PAIR_HPP_

#include "ms/env/match_env.hpp"

namespace toppic {

class EnvPeakPair;
typedef std::shared_ptr<EnvPeakPair> EnvPeakPairPtr;

// Used in dynamic programming algorithm for envelope selection
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
