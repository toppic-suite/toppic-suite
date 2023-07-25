//Copyright (c) 2014 - 2022, The Trustees of Indiana University.
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

#ifndef TOPPIC_TOPFD_ECSCORE_ENV_MS_MAP_ENV_UTIL_HPP
#define TOPPIC_TOPFD_ECSCORE_ENV_MS_MAP_ENV_UTIL_HPP

#include "ms/msmap/ms_map.hpp"
#include "topfd/ecscore/env/seed_env.hpp"
#include "topfd/ecscore/env/ms_map_env.hpp"

namespace toppic {

namespace ms_map_env_util {

double pearsonr(std::vector<double> &X, std::vector<double> &Y);

double compTopThreeInteRatio(SeedEnvPtr seed_ptr, MsMapEnvPtr env_ptr);

double compTopThreeInteRatio(SeedEnvPtr seed_ptr, 
                             std::vector<double> &inte_list);

MsMapEnvPtr getMatchMsMapEnv(MsMapPtr ms_map_ptr, SeedEnvPtr seed_ptr,
                             int sp_id, double mass_tol);

MsMapEnvPtr getMatchMsMapEnv(MsMapPtr ms_map_ptr, SeedEnvPtr seed_ptr,
                             int sp_id, double mass_tole,
                             double min_inte);
}

}

#endif //TOPPIC_UTILITY_FUNCTIONS_HPP
