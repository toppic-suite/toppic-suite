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

std::vector<int> findLocalMinima(std::vector<double> &arr) {
  int n = arr.size();
  std::vector<int> minima;
  for (int i = 1; i < n - 1; i++) {
    if ((arr[i - 1] > arr[i]) and (arr[i] < arr[i + 1])) {
      if (i - 2 > 0)
        if (arr[i - 2] <= arr[i])
          continue;
      if (i + 2 < n)
        if (arr[i + 2] <= arr[i])
          continue;
      minima.push_back(i);
    }
  }
  return minima;
}

std::vector<int> findLocalMaxima(std::vector<double> &arr) {
  int n = arr.size();
  std::vector<int> maxima;
  for (int i = 1; i < n - 1; i++)
    if ((arr[i - 1] < arr[i]) and (arr[i] > arr[i + 1]))
      maxima.push_back(i);
  if (arr[n - 1] > arr[n - 2]) maxima.push_back(n - 1);
  return maxima;
}

SeedEnvPtr getHalfChargeEnv(SeedEnvPtr seed_ptr,
                            double even_odd_peak_ratio) {
  double old_charge = seed_ptr->getCharge();
  if (old_charge < 2) {
    return nullptr;
  }
  int new_charge = int(old_charge / 2);
  double ref_mz = seed_ptr->getReferMz();
  double ref_mass = peak_util::compPeakNeutralMass(ref_mz, new_charge);
  double mono_mass = EnvBase::convertRefMassToMonoMass(ref_mass);
  int sp_id = seed_ptr->getSpecId();
  int env_id = -1;
  double inte = seed_ptr->getSeedInte()/2;
  DeconvPeakPtr peak_ptr = std::make_shared<DeconvPeak>(sp_id, env_id,
                                                        mono_mass, inte,
                                                        new_charge);
  SeedEnvPtr new_seed_ptr = std::make_shared<SeedEnv>(peak_ptr);
  return new_seed_ptr;
}

SeedEnvPtr testHalfChargeState(MsMapPtr matrix_ptr, SeedEnvPtr seed_ptr,
                               EnvSetPtr env_set_ptr, double even_odd_peak_ratio,
                               EcscoreParaPtr para_ptr, double sn_ratio) {
  SeedEnvPtr half_charge_seed = getHalfChargeEnv(seed_ptr, even_odd_peak_ratio);
  bool valid = false;
  valid = seed_env_util::preprocessEnv(matrix_ptr, half_charge_seed, para_ptr, sn_ratio);
  if (!valid)
    return nullptr;
  return half_charge_seed;
}

}
}
