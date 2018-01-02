//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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


#ifndef PROT_FEATURE_VERTEX_A_HPP_
#define PROT_FEATURE_VERTEX_A_HPP_

#include "feature/vertex.hpp"

namespace prot {

class VertexA;
typedef std::shared_ptr<VertexA> VertexAPtr;

class VertexA : public Vertex {
 public:
  VertexA(FeatureMngPtr mng_ptr, int bgn_peak, 
          int pre_win_peak_num, int cur_win_peak_num):
      Vertex(mng_ptr, bgn_peak, pre_win_peak_num, cur_win_peak_num) {
        this_score_ = 0;
        score_ = 0;
        prev_vertex_ = -1;
      }

  VertexA(VertexAPtr ptr):
      Vertex(ptr) {
        this_score_ = ptr->this_score_;
        score_ = ptr->score_;
        prev_vertex_ = ptr->prev_vertex_;
      }

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
