// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef TOPPIC_MS_ENV_MATCH_ENV_REFINE_HPP_
#define TOPPIC_MS_ENV_MATCH_ENV_REFINE_HPP_

#include <vector>

#include "ms/env/match_env.hpp"

namespace toppic {

namespace match_env_refine {

// Re-derive each envelope's theoretical distribution from its reference peak
// (choosing between the current, previous and next monoisotopic position on
// the whole envelope) and rescale the chosen one to the experimental peaks.
// core_ratio is EnvPara::refine_core_ratio_: only theoretical peaks at least
// that fraction of the reference peak, and present in the experimental
// envelope, enter the intensity fit; 0 fits every peak.
void mzRefine(MatchEnvPtrVec& envs, double core_ratio);

void mzRefine(const MatchEnvPtr& env, double core_ratio);

// Indices (into the envelope's peak list) of the core peaks: theoretical
// intensity >= core_ratio * reference intensity and existing in real_env.
// Returns an empty list — meaning "use all peaks" downstream — when fewer
// than 3 such peaks exist, so small envelopes keep the whole-envelope fit.
std::vector<int> coreIdxes(const ExpEnvPtr& real_env, const EnvPtr& theo_env,
                           double core_ratio);

// idxes restricts the fit to those peak indices; empty = all peaks.
void compEnvDist(const EnvPtr& real_env, const EnvPtr& theo_env,
                 const std::vector<int>& idxes, double& best_dist,
                 double& best_ratio);

void compDistWithNorm(const std::vector<double>& real,
                      const std::vector<double>& theo,
                      const std::vector<int>& idxes, double& best_dist,
                      double& best_ratio);

std::vector<double> norm(const std::vector<double>& obs, double ratio);

double compDist(const std::vector<double>& norm,
                const std::vector<double>& theo,
                const std::vector<int>& idxes);

}  // namespace match_env_refine

}  // namespace toppic

#endif
