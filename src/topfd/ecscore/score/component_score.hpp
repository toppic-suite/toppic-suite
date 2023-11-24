//Copyright (c) 2014 - 2023, The Trustees of Indiana University.
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

#ifndef TOPPIC_ECSORE_SCORE_COMPONENT_SCORE_HPP
#define TOPPIC_ECSORE_SCORE_COMPONENT_SCORE_HPP

#include "topfd/ecscore/env_set/env_set.hpp"

namespace toppic {

namespace component_score {

double getAggOddEvenPeakRatio(EnvSetPtr env_set_ptr);

double getAggEnvCorr(EnvSetPtr env_set_ptr);

double get3ScanCorr(EnvSetPtr env_set_ptr, int base_spec, int start_spec);

double getMatchedPeakPercent(EnvSetPtr env_set_ptr,
                             std::vector<std::vector<double>> &theo_map);

int getTheoPeakNum(std::vector<std::vector<double>> &theo_map);

double getConsecutivePeakPercent(EnvSetPtr env_set_ptr);

double getMzErrors(EnvSetPtr env_set_ptr);

}

}

#endif //TOPPIC_GET_COMPONENT_SCORE_HPP
