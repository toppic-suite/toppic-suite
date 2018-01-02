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


#ifndef PROT_FEATURE_DECONV_ONE_SP_HPP_
#define PROT_FEATURE_DECONV_ONE_SP_HPP_

#include "feature/feature_mng.hpp"
#include "feature/deconv_data.hpp"
#include "feature/match_env.hpp"

namespace prot {

class DeconvOneSp {
 public:
  explicit DeconvOneSp(FeatureMngPtr mng_ptr): mng_ptr_(mng_ptr) {}

  void setData(PeakPtrVec &peak_list);

  void setMsLevel(int level) {ms_level_ = level;}

  void setData(PeakPtrVec &peak_list, double max_mass, int max_charge);

  MatchEnvPtrVec getResult() {return result_envs_;}

  void preprocess();

  void run();

  MatchEnvPtrVec postprocess(MatchEnvPtrVec  &dp_envs);

 private:
  FeatureMngPtr mng_ptr_;
  DeconvDataPtr data_ptr_;
  MatchEnvPtrVec result_envs_;
  int ms_level_;
};

typedef std::shared_ptr<DeconvOneSp> DeconvOneSpPtr;

}  // namespace prot
#endif
