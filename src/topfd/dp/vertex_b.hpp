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


#ifndef TOPPIC_TOPFD_DP_VERTEX_B_HPP_
#define TOPPIC_TOPFD_DP_VERTEX_B_HPP_

#include "topfd/dp/vertex.hpp"

namespace toppic {

class VertexB;
typedef std::shared_ptr<VertexB> VertexBPtr;

class VertexB : public Vertex {
 public:
  VertexB(DpParaPtr dp_para_ptr, int bgn_peak, 
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
