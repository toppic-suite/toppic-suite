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

#ifndef TOPPIC_ECSCORE_ENV_COLL_ENV_COLL_UTIL_HPP
#define TOPPIC_ECSCORE_ENV_COLL_ENV_COLL_UTIL_HPP

#include "topfd/ecscore/ecscore_para.hpp"
#include "topfd/ecscore/spectrum/peak_matrix.hpp"
#include "topfd/ecscore/envelope/seed_envelope.hpp"
#include "topfd/ecscore/env_set/env_set.hpp"
#include "topfd/ecscore/env_coll/env_coll.hpp"

namespace toppic {

namespace env_coll_util {

EnvCollPtr findEnvColl(PeakMatrixPtr matrix_ptr, SeedEnvelopePtr seed_ptr,
                       EcscoreParaPtr para_ptr, double sn_ratio);

EnvSetPtrVec getChargeEnvList(PeakMatrixPtr matrix_ptr, SeedEnvelopePtr env_ptr, 
                              EnvSetPtr env_set_ptr, EcscoreParaPtr para_ptr, 
                              double sn_ratio);

bool checkExistingFeatures(PeakMatrixPtr matrix_ptr, EnvCollPtr coll_ptr, 
                           EnvCollPtrVec &env_coll_list, EcscoreParaPtr para_ptr);

}
}


#endif 

