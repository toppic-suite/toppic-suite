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

#ifndef TOPPIC_UTILITY_FUNCTIONS_HPP
#define TOPPIC_UTILITY_FUNCTIONS_HPP

#include <cmath>
#include <vector>
#include <iostream>
#include "topfd/feature_detect/envelope/seed_envelope.hpp"
#include "topfd/feature_detect/env_set/env_set.hpp"
#include "topfd/feature_detect/env_collection/env_collection.hpp"
#include "topfd/feature_detect/envelope/evaluate_envelope.hpp"

namespace toppic {
  namespace utility_functions {
    double pearsonr(std::vector<double> &X, std::vector<double> &Y);

    std::vector<int> findLocalMinima(std::vector<double> &arr);

    std::vector<int> findLocalMaxima(std::vector<double> &arr);

    SeedEnvelope test_half_charge_state(PeakMatrix &peak_matrix, SeedEnvelope &env, EnvSet &top_peak_env_set,
                                        double even_odd_peak_ratios, FeatureParaPtr para_ptr, double sn_ratio);
  }
}

#endif //TOPPIC_UTILITY_FUNCTIONS_HPP
