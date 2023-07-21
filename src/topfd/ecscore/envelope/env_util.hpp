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

#ifndef TOPPIC_TOPFD_ECSCORE_ENVELOPE_ENV_UTIL_HPP
#define TOPPIC_TOPFD_ECSCORE_ENVELOPE_ENV_UTIL_HPP

#include <vector>

#include "ms/msmap/ms_map.hpp"
#include "topfd/ecscore/ecscore_para.hpp"
#include "topfd/ecscore/envelope/seed_env.hpp"
#include "topfd/ecscore/env_set/env_set.hpp"

namespace toppic {

namespace env_util {

double pearsonr(std::vector<double> &X, std::vector<double> &Y);

std::vector<int> findLocalMinima(std::vector<double> &arr);

std::vector<int> findLocalMaxima(std::vector<double> &arr);

SeedEnvPtr testHalfChargeState(MsMapPtr matrix_ptr, SeedEnvPtr seed_ptr,
                               EnvSetPtr env_set_ptr, double even_odd_peak_ratios,
                               EcscoreParaPtr para_ptr, double sn_ratio);

}

}

#endif //TOPPIC_UTILITY_FUNCTIONS_HPP
