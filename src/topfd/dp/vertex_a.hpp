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


#ifndef TOPPIC_TOPFD_DP_VERTEX_A_HPP_
#define TOPPIC_TOPFD_DP_VERTEX_A_HPP_

#include "topfd/dp/vertex.hpp"

namespace toppic {

class VertexA;
typedef std::shared_ptr<VertexA> VertexAPtr;

class VertexA : public Vertex {
 public:
  VertexA(DpParaPtr dp_para_ptr, int bgn_peak, 
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
