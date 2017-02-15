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


#ifndef PROT_FEATURE_VERTEX_A_HPP_
#define PROT_FEATURE_VERTEX_A_HPP_

#include "feature/vertex.hpp"

namespace prot {

class VertexA;
typedef std::shared_ptr<VertexA> VertexAPtr;

class VertexA : public Vertex {
 public:
  VertexA(FeatureMngPtr mng_ptr, int bgn_peak, 
          int pre_win_peak_num, int cur_win_peak_num);
  VertexA(VertexAPtr ptr);

  bool addPreEnv(MatchEnvPtr env, int max_overlap);

  double getThisScr() {return this_score_;}

  double getScrA() {return score_;}

  void setScrA(double score) {score_ = score;}

  int getPreA() {return prev_vertex_;}

  void setPreA(int prev) {prev_vertex_ = prev;}

 private:
  // the sum of score_s of envelopes in previous window 
  double this_score_;
  // current score_ for dp, only the score_s of pre_match_env are included. The
  // reason is that we can only determine the score_ for sharing model when
  // envelopes in the next window are determined.
  double score_;
  // previous vertex for backtracking 
  int prev_vertex_;
};

typedef std::vector<VertexAPtr> VertexAPtrVec;
typedef std::vector<VertexAPtrVec> VertexAPtr2D;

}
#endif
