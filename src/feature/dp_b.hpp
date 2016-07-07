#ifndef PROT_FEATURE_DP_B_HPP_
#define PROT_FEATURE_DP_B_HPP_

#include "feature/dp.hpp"
#include "feature/vertex_b.hpp"

namespace prot {

class DpB : public Dp {
 public:
  DpB(DeconvDataPtr data_ptr, MatchEnvPtr2D &win_envs_, FeatureMngPtr mng_ptr);

  void initGraph();

  void dp();

  void backtrace();

 private:
  VertexBPtr2D vertices_; 
};

}
#endif
