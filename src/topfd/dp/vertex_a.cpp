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


#include "topfd/dp/vertex_a.hpp"

namespace toppic {

VertexA::VertexA(DpParaPtr dp_para_ptr, int bgn_peak, 
                 int pre_win_peak_num, int cur_win_peak_num):
    Vertex(dp_para_ptr, bgn_peak, pre_win_peak_num, cur_win_peak_num) {
      this_score_ = 0;
      score_ = 0;
      prev_vertex_ = -1;
    }

VertexA::VertexA(VertexAPtr ptr):
    Vertex(ptr) {
      this_score_ = ptr->this_score_;
      score_ = ptr->score_;
      prev_vertex_ = ptr->prev_vertex_;
    }

bool VertexA::addPreEnv(MatchEnvPtr env, int max_overlap) {
  bool result = Vertex::addPreEnv(env, max_overlap);
  this_score_ += env->getScore();
  return result;
}

}
