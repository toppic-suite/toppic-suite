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


#include "topfd/dp/dp.hpp"

namespace toppic {

Dp::Dp (DeconvDataPtr data_ptr, MatchEnvPtr2D &win_envs, 
        DpParaPtr dp_para_ptr, double score_error_tolerance):
    data_ptr_(data_ptr),
    dp_para_ptr_(dp_para_ptr), 
    score_error_tolerance_(score_error_tolerance),
    win_envs_(win_envs) {
      win_num_ = data_ptr_->getWinNum();
    }

// add an envelope list to result list 
void Dp::addEnv(MatchEnvPtrVec &result, MatchEnvPtrVec &prev_env) {
  result.insert(result.end(), prev_env.begin(), prev_env.end());
}

}
