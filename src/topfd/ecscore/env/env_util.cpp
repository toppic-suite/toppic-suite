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

#include "ms/spec/peak_util.hpp"
#include "topfd/ecscore/env/seed_env_util.hpp"
#include "topfd/ecscore/env/env_util.hpp"

namespace toppic {

namespace env_util {

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

}

}
