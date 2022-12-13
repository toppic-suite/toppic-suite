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

#ifndef TOPPIC_GET_COMPONENT_SCORE_HPP
#define TOPPIC_GET_COMPONENT_SCORE_HPP

#include <valarray>
#include "env_util.hpp"
#include "topfd/envcnn/env_cnn.hpp"
#include "topfd/feature_detect/env_collection/env_collection.hpp"

namespace toppic {
  namespace component_score {
    double get_agg_odd_even_peak_ratio(EnvSet &env_set);

    double get_agg_env_corr(EnvSet &env_set);

    double get_3_scan_corr(EnvSet &env_set, int base_spec, int start_spec);

    double get_matched_peaks_percent(EnvSet &env_set, std::vector<std::vector<double>> theo_map);

    int get_num_theo_peaks(std::vector<std::vector<double>> &theo_map);

    double get_consecutive_peaks_percent(EnvSet &env_set);

    double get_mz_errors(EnvSet &env_set);
  }
}

#endif //TOPPIC_GET_COMPONENT_SCORE_HPP
