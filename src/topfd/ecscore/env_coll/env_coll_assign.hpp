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

#ifndef TOPPIC_TOPFD_ECSCORE_ENV_COLL_ENV_COLL_ASSIGN_HPP
#define TOPPIC_TOPFD_ECSCORE_ENV_COLL_ENV_COLL_ASSIGN_HPP

#include "ms/feature/frac_feature.hpp"
#include "ms/feature/spec_feature.hpp"
#include "ms/msmap/ms_map.hpp"

#include "topfd/common/topfd_para.hpp"
#include "topfd/ecscore/para/ecscore_para.hpp"
#include "topfd/ecscore/score/ecscore.hpp"

#include "topfd/ecscore/env_coll/env_coll.hpp"

namespace toppic {

namespace env_coll_assign {

void assignEnvColls(FracFeaturePtrVec &frac_feature_list,
                    EnvCollPtrVec &env_coll_list,
                    MsHeaderPtr2D &ms2_header_ptr_2d,
                    SpecFeaturePtrVec &ms2_feature_list, 
                    double score_cutoff); 

bool checkEnvColl(MsHeaderPtr header_ptr, 
                  EnvCollPtrVec &env_coll_list); 


}

}

#endif 
