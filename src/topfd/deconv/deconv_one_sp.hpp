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

#ifndef TOPPIC_TOPFD_DECONV_ONE_SP_HPP_
#define TOPPIC_TOPFD_DECONV_ONE_SP_HPP_

#include "topfd/spec/deconv_data.hpp"
#include "ms/env/env_para.hpp"
#include "ms/env/match_env.hpp"
#include "topfd/dp/dp_para.hpp"

namespace toppic {

class DeconvOneSp {
 public:
  explicit DeconvOneSp(EnvParaPtr env_para_ptr, DpParaPtr dp_para_ptr): 
      env_para_ptr_(env_para_ptr), dp_para_ptr_(dp_para_ptr) {}

  void setData(PeakPtrVec &peak_list);

  void setMsLevel(int level) {ms_level_ = level;}

  void setData(PeakPtrVec &peak_list, double max_mass, int max_charge);

  MatchEnvPtrVec getResult() {return result_envs_;}

  void preprocess();

  void run();

  MatchEnvPtrVec postprocess(MatchEnvPtrVec  &dp_envs);

  EnvParaPtr getEnvPara(){return env_para_ptr_;}
  DpParaPtr getDpPara(){return dp_para_ptr_;}

 private:
  EnvParaPtr env_para_ptr_;
  DpParaPtr dp_para_ptr_;
  DeconvDataPtr data_ptr_;
  MatchEnvPtrVec result_envs_;
  int ms_level_;
};

typedef std::shared_ptr<DeconvOneSp> DeconvOneSpPtr;

}  // namespace toppic
#endif
