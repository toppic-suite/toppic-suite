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


#ifndef PROT_FEATURE_DP_HPP_
#define PROT_FEATURE_DP_HPP_

#include "feature/deconv_data.hpp"
#include "feature/feature_mng.hpp"
#include "feature/match_env.hpp"

namespace toppic {

class Dp {
 public:
  Dp(DeconvDataPtr data_ptr, MatchEnvPtr2D &win_envs, 
     FeatureMngPtr mng_ptr);
  void addEnv(MatchEnvPtrVec &result, MatchEnvPtrVec &prev_env);

  MatchEnvPtrVec getResult() {return results_;}

 protected:
  DeconvDataPtr data_ptr_;
  FeatureMngPtr mng_ptr_;
  MatchEnvPtr2D win_envs_;
  int win_num_;
  MatchEnvPtrVec results_;
};

}
#endif
