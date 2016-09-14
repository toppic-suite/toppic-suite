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


#ifndef PROT_FEATURE_VERTEX_HPP_
#define PROT_FEATURE_VERTEX_HPP_

#include "feature/feature_mng.hpp"
#include "feature/match_env.hpp"
#include "feature/env_peak_pair.hpp"

namespace prot {

class Vertex;
typedef std::shared_ptr<Vertex> VertexPtr;

class Vertex {
 public:
  Vertex(FeatureMngPtr mng_ptr, int bgn_peak, int pre_win_peak_num, int cur_win_peak_num);

  Vertex(VertexPtr ptr);

  Vertex() {};

  virtual bool addPreEnv(MatchEnvPtr env, int max_overlap);
  bool addCurEnv(MatchEnvPtr env, int max_overlap);
  void trim();
  static bool checkConsist(VertexPtr pre, VertexPtr cur, int max_env_per_peak);

  int getMatchEnvSize() {return prec_match_envs_.size() + cur_match_envs_.size();}

  std::vector<int> getPeakUseCnts() {return peak_use_cnts_;}

  MatchEnvPtrVec getPreMatchEnvs() {return prec_match_envs_;}

  MatchEnvPtrVec getCurMatchEnvs() {return cur_match_envs_;}

  int getPreMatchEnvSize() {return prec_match_envs_.size();}

  static double getShareScr(VertexPtr pre, VertexPtr cur, double error_tole);

 protected:
  FeatureMngPtr mng_ptr_;
  int bgn_peak_;
  int peak_num_;
  int cur_bgn_pos_;

  std::vector<int> peak_use_cnts_;

  // for each envelope, there is a peak list 
  MatchEnvPtrVec prec_match_envs_;
  MatchEnvPtrVec cur_match_envs_;

  // There is a list of EnvPeakPairs for each peak 
  EnvPeakPairPtr2D pre_env_peak_pairs_;
  EnvPeakPairPtr2D cur_env_peak_pairs_;

  bool passDblIncrCheck(MatchEnvPtr env);

  bool passDblIncrCheck(MatchEnvPtrVec &env_list); 
};

typedef std::vector<VertexPtr> VertexPtrVec;

}
#endif
