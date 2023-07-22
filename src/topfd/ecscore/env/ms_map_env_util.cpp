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

#include <cmath>

#include "common/util/logger.hpp"

#include "topfd/ecscore/env/ms_map_env_util.hpp"

namespace toppic {

namespace ms_map_env_util {

double pearsonr(std::vector<double> &X, std::vector<double> &Y) {
  int n = X.size();
  double sum_X = 0, sum_Y = 0, sum_XY = 0;
  double squareSum_X = 0, squareSum_Y = 0;
  for (int i = 0; i < n; i++) {
    sum_X = sum_X + X[i];
    sum_Y = sum_Y + Y[i];
    sum_XY = sum_XY + X[i] * Y[i];
    squareSum_X = squareSum_X + X[i] * X[i];
    squareSum_Y = squareSum_Y + Y[i] * Y[i];
  }
  double corr = (double) (n * sum_XY - sum_X * sum_Y) /
    sqrt((n * squareSum_X - sum_X * sum_X) * (n * squareSum_Y - sum_Y * sum_Y));
  return corr;
}

double compTopThreeInteRatio(SeedEnvPtr seed_ptr, MsMapEnvPtr env_ptr) {
  double seed_inte = seed_ptr->compTopThreeInteSum();
  int ref_idx = seed_ptr->getReferIdx();
  double env_inte = env_ptr->compTopThreeInteSum(ref_idx);
  if (seed_inte == 0) {
    LOG_WARN("Empty peak list!");
    return 0;
  }
  return env_inte/seed_inte;
}

}

}
