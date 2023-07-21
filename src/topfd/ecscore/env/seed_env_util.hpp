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

#ifndef TOPPIC_TOPFD_ECSCORE_ENVELOPE_SEED_ENV_UTIL_HPP
#define TOPPIC_TOPFD_ECSCORE_ENVELOPE_SEED_ENV_UTIL_HPP

#include "ms/msmap/ms_map.hpp"
#include "topfd/ecscore/ecscore_para.hpp"
#include "topfd/ecscore/env/seed_env.hpp"

namespace toppic {

namespace seed_env_util {

bool preprocessEnv(MsMapPtr matrix_ptr, SeedEnvPtr seed_ptr,
                   EcscoreParaPtr para_ptr, double sn_ratio);

bool evalEnv(MsMapPtr matrix_ptr, SeedEnvPtr seed_ptr,
             double mass_tol, double corr_tol, double snr);

bool simplePreprocessEnv(MsMapPtr matrix_ptr, SeedEnvPtr seed_ptr,
                         EcscoreParaPtr para_ptr, double sn_ratio); 

bool testChargeState(int charge, std::vector<double> &seed_env_inte);

bool evalEnvPair(std::vector<double> &exp_env_inte, std::vector<double> &theo_inte, double corr_tol);

}

}


#endif //TOPPIC_EVALUATE_ENVELOPE_HPP
