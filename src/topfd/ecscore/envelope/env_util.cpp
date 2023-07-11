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

#include "ms/spec/peak_util.hpp"
#include "topfd/ecscore/envelope/seed_env_util.hpp"
#include "topfd/ecscore/envelope/env_util.hpp"

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

SeedEnvelopePtr getHalfChargeEnv(SeedEnvelopePtr seed_ptr, 
                                 double even_odd_peak_ratio) {
  double mass = seed_ptr->getMass();
  int charge = seed_ptr->getCharge();
  double mz = peak_util::compMz(mass, charge);
  std::vector<double> distribution = seed_ptr->getPosList();
  if (even_odd_peak_ratio < 0)
    mz = mz + (distribution[1] - distribution[0]);
  int new_charge = int(charge / 2);
  if (new_charge == 0)
    new_charge = new_charge + 1;
  double new_mass = peak_util::compPeakNeutralMass(mz, new_charge);
  // get a reference distribution based on the base mass
  EnvPtr ref_env_ptr = EnvBase::getEnvByMonoMass(new_mass);
  if (ref_env_ptr == nullptr) return nullptr; 
  EnvPtr theo_env_ptr = ref_env_ptr->distrToTheoMono(mz, new_charge);
  std::vector<double> env_peaks_mz, env_peaks_inte;
  for (int i = 0; i < theo_env_ptr->getPeakNum(); i++) {
    env_peaks_mz.push_back(theo_env_ptr->getMz(i));
    env_peaks_inte.push_back(theo_env_ptr->getIntensity(i));
  }
  SeedEnvelopePtr sp_peak =
    std::make_shared<SeedEnvelope>(seed_ptr->getSpecId(), seed_ptr->getEnvId(), 
                                   seed_ptr->getPos(), new_mass, seed_ptr->getInte(),
                                      new_charge, env_peaks_mz, env_peaks_inte);
  return sp_peak;
}

SeedEnvelopePtr testHalfChargeState(PeakMatrixPtr matrix_ptr, SeedEnvelopePtr seed_ptr, 
                                    EnvSetPtr env_set_ptr, double even_odd_peak_ratio, 
                                    EcscoreParaPtr para_ptr, double sn_ratio) {
  SeedEnvelopePtr half_charge_seed = getHalfChargeEnv(seed_ptr, even_odd_peak_ratio);
  bool valid = false;
  valid = seed_env_util::preprocessEnv(matrix_ptr, half_charge_seed, para_ptr, sn_ratio);
  if (!valid)
    return nullptr;
  return half_charge_seed;
}

}
}
