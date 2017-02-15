// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#ifndef PROT_FEATURE_MATCH_ENV_UTIL_HPP_
#define PROT_FEATURE_MATCH_ENV_UTIL_HPP_

#include "spec/peak.hpp"
#include "feature/match_env.hpp"

namespace prot {

class MatchEnvUtil {
 public:
  static MatchEnvPtrVec sortOnMz(MatchEnvPtrVec &ori_envs);

  static MatchEnvPtrVec sortOnMass(MatchEnvPtrVec &ori_envs);

  static std::vector<double> getMassList(MatchEnvPtrVec &envs);

  static std::vector<int> getChargeList(MatchEnvPtrVec &envs);

  static std::vector<double> getChargeOneMassList(MatchEnvPtrVec &envs);

  static std::vector<double> getIntensitySums(MatchEnvPtrVec &envs);

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
