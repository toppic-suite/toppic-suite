#include "feature/vertex_a.hpp"

namespace prot {

VertexA::VertexA(FeatureMngPtr mng_ptr, int bgn_peak, 
                 int pre_win_peak_num, int cur_win_peak_num):
    Vertex(mng_ptr, bgn_peak, pre_win_peak_num, cur_win_peak_num) {
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
