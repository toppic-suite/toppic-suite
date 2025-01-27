//Copyright (c) 2014 - 2025, The Trustees of Indiana University, Tulane University.
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

#ifndef TOPPIC_TOPFD_ECSCORE_ENV_SEED_ENV_UTIL_HPP
#define TOPPIC_TOPFD_ECSCORE_ENV_SEED_ENV_UTIL_HPP

#include "ms/msmap/ms_map.hpp"
#include "topfd/ecscore/para/ecscore_para.hpp"
#include "topfd/ecscore/env/seed_env.hpp"

namespace toppic {

namespace seed_env_util {

double pearsonr(std::vector<double> &X, std::vector<double> &Y); 

SeedEnvPtr preprocessSeedEnvPtr(SeedEnvPtr seed_ptr, MsMapPtr ms_map_ptr,
                                EcscoreParaPtr para_ptr, double sn_ratio);

SeedEnvPtr relaxProcessSeedEnvPtr(SeedEnvPtr seed_ptr, MsMapPtr ms_map_ptr,
                                  EcscoreParaPtr para_ptr, double sn_ratio);

SeedEnvPtr testHalfChargeEnv(SeedEnvPtr seed_ptr, MsMapPtr ms_map_ptr, 
                             double even_odd_log_ratio, EcscoreParaPtr para_ptr, 
                             double sn_ratio); 

}

}

#endif 
