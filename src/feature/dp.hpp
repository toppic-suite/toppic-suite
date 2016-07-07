#ifndef PROT_FEATURE_DP_HPP_
#define PROT_FEATURE_DP_HPP_

#include "feature/deconv_data.hpp"
#include "feature/feature_mng.hpp"
#include "feature/match_env.hpp"

namespace prot {

class Dp {
 public:
  Dp(DeconvDataPtr data_ptr, MatchEnvPtr2D &win_envs, 
     FeatureMngPtr mng_ptr);
  void addEnv(MatchEnvPtrVec &result, MatchEnvPtrVec &prev_env);

  MatchEnvPtrVec getResult() {return results_;}

 private:
  DeconvDataPtr data_ptr_;
  FeatureMngPtr mng_ptr_;
  MatchEnvPtr2D win_envs_;
  int win_num_;
  MatchEnvPtrVec results_;
};

}
#endif
