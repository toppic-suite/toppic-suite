// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
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


#ifndef PROT_FEATURE_MATCH_ENVELOPE_HPP_
#define PROT_FEATURE_MATCH_ENVELOPE_HPP_

#include <memory>
#include <vector>

#include "feature/feature_mng.hpp" 
#include "feature/envelope.hpp" 
#include "feature/real_env.hpp" 

namespace prot {

class MatchEnv;

typedef std::shared_ptr<MatchEnv> MatchEnvPtr;

class MatchEnv {
 public:
  MatchEnv(int mass_group, EnvelopePtr theo_env_ptr, RealEnvPtr real_env_ptr);

  void compScr(FeatureMngPtr mng_ptr);

  static bool cmpScoreDec(const MatchEnvPtr &a, const MatchEnvPtr &b) { 
    return a->getScore() > b->getScore();
  }

  double calcPeakScr(int id_x, double inte_sum, double tolerance);

  int getId() {return id_;}

  int getMassGroup() {return mass_group_;}

  RealEnvPtr getRealEnvPtr() {return real_env_ptr_;}

  EnvelopePtr getTheoEnvPtr() {return theo_env_ptr_;}

  double getScore() {return score_;}

  void setId(int id) {id_ = id;}

  void setTheoEnvPtr(EnvelopePtr theo_env_ptr) {theo_env_ptr_ = theo_env_ptr;}

 private:
  int id_;
  // we divide envelopes into several groups based on monoisotopic  masses  
  int mass_group_;
  double score_;
  EnvelopePtr theo_env_ptr_;
  RealEnvPtr real_env_ptr_;

  double calcShareInteAccu(int id_x, double inte_sum);

  double calcMzFactor(int id_x, double shift, double tolerance);

  double calcIntensityFactor(double theo_inte, double real_inte);

  double calcIntensityFactor(int id_x, double ratio);

  double findBestShift(FeatureMngPtr mng_ptr);

  double findBestRatio(FeatureMngPtr mng_ptr);

  double calcScrWithSftRatio(double shift, double ratio, double tolerance);
};

typedef std::vector<MatchEnvPtr> MatchEnvPtrVec;
typedef std::vector<MatchEnvPtrVec> MatchEnvPtr2D;

}

#endif
