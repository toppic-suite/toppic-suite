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

#ifndef TOPPIC_TOPFD_DECONV_DECONV_SINGLE_SP_HPP_
#define TOPPIC_TOPFD_DECONV_DECONV_SINGLE_SP_HPP_

#include "ms/env/env_para.hpp"
#include "topfd/dp/dp_para.hpp"
#include "ms/env/match_env.hpp"
#include "topfd/spec/deconv_data.hpp"

namespace toppic {

class DeconvSingleSp {
 public:
  explicit DeconvSingleSp(TopfdParaPtr topfd_para_ptr, PeakPtrVec &peak_list,
                          int ms_level, double max_mass, int max_charge);
  MatchEnvPtrVec deconv();
  
  void postprocess(MatchEnvPtrVec  &dp_envs);

 private:
  TopfdParaPtr topfd_para_ptr_;
  EnvParaPtr env_para_ptr_;
  DpParaPtr dp_para_ptr_;
  // spectral data
  DeconvDataPtr data_ptr_;
  // resulting envelopes
  MatchEnvPtrVec result_envs_;
};

typedef std::shared_ptr<DeconvSingleSp> DeconvSingleSpPtr;

}

#endif

