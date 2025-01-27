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

#ifndef TOPPIC_TOPFD_ECSCORE_SCORE_COMP_ENV_CNN_SCORE_HPP
#define TOPPIC_TOPFD_ECSCORE_SCORE_COMP_ENV_CNN_SCORE_HPP

#include "ms/msmap/ms_map.hpp"
#include "topfd/ecscore/env_coll/env_coll.hpp"

namespace toppic {

namespace comp_env_cnn_score {

double compEnvcnnScore(MsMapPtr matrix_ptr, EnvCollPtr env_coll_ptr);

}
}

#endif 
