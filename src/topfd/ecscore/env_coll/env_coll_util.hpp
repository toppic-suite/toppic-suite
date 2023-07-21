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

#include "ms/spec/deconv_ms.hpp"
#include "ms/feature/frac_feature.hpp"

#include "ms/msmap/ms_map.hpp"
#include "topfd/ecscore/ecscore_para.hpp"
#include "topfd/ecscore/env/seed_env.hpp"
#include "topfd/ecscore/env_set/env_set.hpp"
#include "topfd/ecscore/env_coll/env_coll.hpp"

namespace toppic {

namespace env_coll_util {

EnvCollPtr findEnvColl(MsMapPtr matrix_ptr, SeedEnvPtr seed_ptr,
                       EcscoreParaPtr para_ptr, double sn_ratio);

EnvCollPtr findEnvCollWithSingleEnv(MsMapPtr matrix_ptr, SeedEnvPtr seed_ptr,
                                    EcscoreParaPtr para_ptr, double sn_ratio); 

EnvSetPtrVec getChargeEnvList(MsMapPtr matrix_ptr, SeedEnvPtr seed_ptr,
                              EnvSetPtr seed_env_set_ptr, EcscoreParaPtr para_ptr,
                              double sn_ratio);

bool checkExistingFeatures(MsMapPtr matrix_ptr, EnvCollPtr coll_ptr,
                           EnvCollPtrVec &env_coll_list, EcscoreParaPtr para_ptr);

FracFeaturePtr getFracFeature(int feat_id, DeconvMsPtrVec &ms1_ptr_vec, int
                              frac_id, std::string &file_name,
                              EnvCollPtr coll_ptr, MsMapPtr matrix_ptr, double snr);

}
}


#endif 

