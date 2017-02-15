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


#ifndef PROT_FEATURE_VERTEX_B_HPP_
#define PROT_FEATURE_VERTEX_B_HPP_

#include "feature/vertex.hpp"

namespace prot {

class VertexB;
typedef std::shared_ptr<VertexB> VertexBPtr;

class VertexB : public Vertex {
 public:
  VertexB(FeatureMngPtr mng_ptr, int bgn_peak, 
          int pre_win_peak_num, int cur_win_peak_num, int env_num);
  VertexB(VertexBPtr ptr);

  bool addPreEnv(MatchEnvPtr env, int max_overlap) {return false;}

  double getScoreB(int i) {return scores_[i];}

  int getPrevB(int i) {return prevs_[i];}

  void setScoreB(int i, double s) {scores_[i] = s;}

  void setPrevB(int i, int p) {prevs_[i] = p;}

 private:
    // dp scores for different number of envelopes 
  std::vector<double> scores_;
    // previous vertex id 
  std::vector<int> prevs_;
};

typedef std::vector<VertexBPtr> VertexBPtrVec;
typedef std::vector<VertexBPtrVec> VertexBPtr2D;

}
#endif
