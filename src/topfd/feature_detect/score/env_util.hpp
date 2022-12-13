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

#ifndef TOPPIC_ENV_UTIL_HPP
#define TOPPIC_ENV_UTIL_HPP

#include <iostream>
#include <vector>
#include "topfd/feature_detect/env_set/env_set.hpp"

namespace toppic {
  namespace env_utils {
    double get_mz(double mass, int charge);

    double get_mass(double mz, int charge);

    std::vector<double> get_aggregate_envelopes_mz(EnvSet &env_set);

    std::vector<double> get_aggregate_envelopes_inte(EnvSet &env_set);

    double calcInteRatio_scan(std::vector<double> &theo_envelope_inte, std::vector<double> &exp_envelope_inte);
  }
}
#endif //TOPPIC_ENV_UTIL_HPP
