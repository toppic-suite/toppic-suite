//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#ifndef PROT_FEATURE_ENV_DETECT_HPP_
#define PROT_FEATURE_ENV_DETECT_HPP_

#include <memory>
#include <vector>

#include "feature/envelope.hpp"
#include "feature/match_env.hpp"
#include "feature/feature_mng.hpp"
#include "feature/deconv_data.hpp"

namespace prot {

class EnvDetect {
 public:
  static double calcInteRatio(EnvelopePtr theo_env, PeakPtrVec &peak_list, 
                              double tolerance);
  static MatchEnvPtr detectEnv(PeakPtrVec &peak_list, int base_peak,
                               int charge, double max_mass, FeatureMngPtr mng_ptr);

  static MatchEnvPtr detectEnv(PeakPtrVec &peak_list, double mono_mass, 
                               int charge, FeatureMngPtr mng_ptr);

  static MatchEnvPtr2D getCandidate(DeconvDataPtr data_ptr, FeatureMngPtr mng_ptr);
};

}

#endif
