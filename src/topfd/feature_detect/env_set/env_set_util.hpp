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

#ifndef TOPPIC_ENV_SET_UTIL_HPP
#define TOPPIC_ENV_SET_UTIL_HPP

#include <valarray>
#include "exp_envelope.hpp"
#include "topfd/feature_detect/env_set/env_set.hpp"
#include "topfd/feature_detect/spectrum/peak_row.hpp"
#include "topfd/feature_detect/envelope/simple_peak.hpp"
#include "topfd/feature_detect/envelope/seed_envelope.hpp"
#include "topfd/feature_detect/envelope/evaluate_envelope.hpp"
#include "topfd/feature_detect/spectrum/peak_matrix.hpp"

namespace toppic {
namespace env_set_util {
  EnvSet get_env_set(PeakMatrix& peak_matrix, SeedEnvelope& seed_env, FeatureParaPtr para_ptr, double sn_ratio);
  EnvSet find_env_set(PeakMatrix& peak_matrix, SeedEnvelope& seed_env, int start_spec_id, int end_spec_id, FeatureParaPtr para_ptr, double sn_ratio);
  void comp_peak_start_end_idx(PeakMatrix& peak_matrix, SeedEnvelope &seed_env, double error_tole);
  void remove_non_match_envs(std::vector<ExpEnvelope>& env_list, int refer_idx);
  bool check_valid_env_set(PeakMatrix& peak_matrix, EnvSet& env_set);
  bool check_valid_env_set_seed_env(PeakMatrix& peak_matrix, EnvSet& env_set, int max_miss_peak);
  bool check_valid_env_set_seed_env_sparse(PeakMatrix& peak_matrix, EnvSet& env_set, int max_miss_peak);
  ExpEnvelope get_match_exp_env(PeakMatrix& peak_matrix, SeedEnvelope& seed_env, int sp_id, double mass_tol);
}
}


#endif //TOPPIC_ENV_SET_UTIL_HPP
