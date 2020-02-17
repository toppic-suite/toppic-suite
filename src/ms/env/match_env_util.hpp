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


#ifndef TOPPIC_TOPFD_ENV_MATCH_ENV_UTIL_HPP_
#define TOPPIC_TOPFD_ENV_MATCH_ENV_UTIL_HPP_

#include "ms/spec/peak.hpp"
#include "ms/spec/ms_header.hpp"
#include "ms/spec/deconv_ms.hpp"
#include "ms/env/match_env.hpp"

namespace toppic {

namespace match_env_util {

std::vector<double> getMassList(const MatchEnvPtrVec &envs);

std::vector<int> getChargeList(const MatchEnvPtrVec &envs);

std::vector<double> getChargeOneMassList(const MatchEnvPtrVec &envs);

std::vector<double> getIntensitySums(const MatchEnvPtrVec &envs);

void assignIntensity(PeakPtrVec &ms, MatchEnvPtrVec &envs);

PeakPtrVec rmAnnoPeak(PeakPtrVec &ms, MatchEnvPtrVec &envs);

MatchEnvPtrVec addLowMassPeak(MatchEnvPtrVec &envs, std::vector<PeakPtr> &ms, double tolerance);

MatchEnvPtr getNewMatchEnv(PeakPtrVec &ms, int idx, double tolerance);

MatchEnvPtrVec addMultipleMass(MatchEnvPtrVec &envs, MatchEnvPtr2D &candidates,
                               double multi_min_mass, int multi_min_charge, double min_ratio);

DeconvMsPtr getDeconvMsPtr(MsHeaderPtr header_ptr, MatchEnvPtrVec &envs);

}  // namespace match_env_util

}  // namespace toppic

#endif
