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


#ifndef PROT_FEATURE_MATCH_ENV_UTIL_HPP_
#define PROT_FEATURE_MATCH_ENV_UTIL_HPP_

#include "spec/peak.hpp"
#include "feature/match_env.hpp"

namespace prot {

class MatchEnvUtil {
 public:
  static std::vector<double> getMassList(const MatchEnvPtrVec &envs);

  static std::vector<int> getChargeList(const MatchEnvPtrVec &envs);

  static std::vector<double> getChargeOneMassList(const MatchEnvPtrVec &envs);

  static std::vector<double> getIntensitySums(const MatchEnvPtrVec &envs);

  static void assignIntensity(PeakPtrVec &ms, MatchEnvPtrVec &envs);

  static PeakPtrVec rmAnnoPeak(PeakPtrVec &ms, MatchEnvPtrVec &envs);

  static MatchEnvPtrVec addLowMassPeak(MatchEnvPtrVec &envs, std::vector<PeakPtr> &ms, 
                                       double tolerance);

  static MatchEnvPtr getNewMatchEnv(PeakPtrVec &ms, int idx, double tolerance);

  static MatchEnvPtrVec addMultipleMass(MatchEnvPtrVec &envs, MatchEnvPtr2D &candidates,
                                        double multi_min_mass, int multi_min_charge, double min_ratio);

};

}

#endif
