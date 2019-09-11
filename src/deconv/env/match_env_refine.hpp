//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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


#ifndef TOPPIC_DECONV_ENV_MATCH_ENV_REFINE_HPP_
#define TOPPIC_DECONV_ENV_MATCH_ENV_REFINE_HPP_

#include "deconv/env/match_env.hpp"

namespace toppic {

namespace match_env_refine {

void mzRefine(MatchEnvPtrVec &envs);

void mzRefine(MatchEnvPtr env);

void compEnvDist(EnvelopePtr real_env, EnvelopePtr theo_env,
                 double &best_dist, double &best_ratio);

void compDistWithNorm(const std::vector<double>& real,
                      const std::vector<double>& theo, 
                      double &best_dist, double &best_ratio);

std::vector<double> norm(const std::vector<double> &obs, double ratio);

double compDist(const std::vector<double> &norm, const std::vector<double> &theo);

} // match_env_refine

}  // namespace toppic

#endif
