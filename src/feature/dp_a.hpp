#ifndef PROT_FEATURE_DP_A_HPP_
#define PROT_FEATURE_DP_A_HPP_

#include "feature/dp.hpp"
#include "feature/vertex_a.hpp"

namespace prot {

class DpA : public Dp {
 public:
  DpA(DeconvDataPtr data_ptr, MatchEnvPtr2D &win_envs_, FeatureMngPtr mng_ptr);

  void initGraph();

  void dp();

  void backtrace();

 private:
  VertexAPtr2D vertices_; 
};

}
#endif
