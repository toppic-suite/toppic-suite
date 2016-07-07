#include "feature/dp.hpp"

namespace prot {

Dp::Dp (DeconvDataPtr data_ptr, MatchEnvPtr2D &win_envs, 
        FeatureMngPtr mng_ptr):
    data_ptr_(data_ptr),
    mng_ptr_(mng_ptr), 
    win_envs_(win_envs) {
      win_num_ = data_ptr_->getWinNum();
    }

// add an envelope list to result list 
void Dp::addEnv(MatchEnvPtrVec &result, MatchEnvPtrVec &prev_env) {
  result.insert(result.end(), prev_env.begin(), prev_env.end());
}

}
