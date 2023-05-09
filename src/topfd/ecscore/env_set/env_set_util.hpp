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

#ifndef TOPPIC_TOPFD_ECSCORE_ENV_SET_ENV_SET_UTIL_HPP
#define TOPPIC_TOPFD_ECSCORE_ENV_SET_ENV_SET_UTIL_HPP

#include <vector>

#include "topfd/ecscore/ecscore_para.hpp"
#include "topfd/ecscore/spectrum/peak_matrix.hpp"
#include "topfd/ecscore/env_set/env_set.hpp"

namespace toppic {

namespace env_set_util {

std::vector<double> getAggregateEnvelopeMz(EnvSetPtr env_set_ptr);

std::vector<double> getAggregateEnvelopeInte(EnvSetPtr env_set_ptr);

double calcInteRatio(std::vector<double> &theo_envelope_inte, std::vector<double> &exp_envelope_inte);

EnvSetPtr getEnvSet(PeakMatrixPtr matrix_ptr, SeedEnvelopePtr seed_ptr, 
                    EcscoreParaPtr para_ptr, double sn_ratio);

EnvSetPtr findEnvSet(PeakMatrixPtr matrix_ptr, SeedEnvelopePtr seed_ptr, 
                     int start_spec_id, int end_spec_id, 
                     EcscoreParaPtr para_ptr, double sn_ratio);

void compPeakStartEndIdx(PeakMatrixPtr matrix_ptr, SeedEnvelopePtr seed_ptr, double error_tole);

bool checkValidEnvSet(PeakMatrixPtr matrix_ptr, EnvSetPtr env_set_ptr);

bool checkValidEnvSetSeedEnv(PeakMatrixPtr matrix_ptr, EnvSetPtr env_set_ptr, int max_miss_peak);

bool checkValidEnvSetSeedEnvSparse(PeakMatrixPtr matrix_ptr, EnvSetPtr env_set_ptr, int max_miss_peak);

ExpEnvelopePtr getMatchExpEnv(PeakMatrixPtr matrix_ptr, SeedEnvelopePtr seed_ptr, int sp_id, double mass_tol);

}

}
#endif 
