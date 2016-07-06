#include "feature/vertex_b.hpp"

namespace prot {

VertexB::VertexB(FeatureMngPtr mng_ptr, int bgn_peak, 
                 int pre_win_peak_num, int cur_win_peak_num, int env_num):
    Vertex(mng_ptr, bgn_peak, pre_win_peak_num, cur_win_peak_num) {
      // zero envelope is considered, so we have env_num + 1 
      scores_.resize(env_num + 1, 0);
      prevs_.resize(env_num + 1, -1);
    }

VertexB::VertexB(VertexBPtr ptr):
    Vertex(ptr) {
    scores_ = ptr->scores_;
    prevs_ = ptr->prevs_;
    }

}
