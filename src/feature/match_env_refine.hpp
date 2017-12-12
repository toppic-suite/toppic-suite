//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#ifndef PROT_FEATURE_MATCH_ENV_REFINE_HPP_
#define PROT_FEATURE_MATCH_ENV_REFINE_HPP_

#include <vector>

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

  static double compDistWithNorm(const std::vector<double>& real,
                                 const std::vector<double>& theo);

  static std::vector<double> norm(const std::vector<double> &obs, double ratio);

  static double compDist(const std::vector<double> &norm, const std::vector<double> &theo);
};

}  // namespace prot

#endif
