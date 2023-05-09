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

#ifndef TOPPIC_ENV_COLL_UTIL_HPP
#define TOPPIC_ENV_COLL_UTIL_HPP

#include "env_collection.hpp"
#include "topfd/feature_detect/env_set/env_set_util.hpp"
#include "topfd/feature_detect/test_output_functions/write_out_files.hpp"
#include "topfd/feature_detect/score/get_component_score.hpp"
#include "topfd/feature_detect/util/utility_functions.hpp"
#include "topfd/feature_detect/spectrum/peak_matrix.hpp"
#include "topfd/feature_detect/envelope/seed_envelope.hpp"
#include "ms/env/env_base.hpp"

namespace toppic {
  namespace env_coll_util {
    EnvCollection find_env_collection(PeakMatrix &peak_matrix, SeedEnvelope &seed_env, FeatureParaPtr para_ptr, double sn_ratio);
    std::vector<EnvSet> get_charge_env_list(PeakMatrix &peak_matrix, SeedEnvelope &env, EnvSet &top_peak_env_set, FeatureParaPtr para_ptr, double sn_ratio);
    bool check_in_existing_features(PeakMatrix &peakMatrix, EnvCollection &env_coll, std::vector<EnvCollection> &env_coll_list, FeatureParaPtr para_ptr);
  }
}


#endif //TOPPIC_ENV_COLL_UTIL_HPP
