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

#include "topfd/ecscore/env_set/env_set_util.hpp"

namespace toppic {

namespace env_set_util {

std::vector<double> getAggregateEnvelopeInte(EnvSetPtr env_set_ptr) {
  ExpEnvelopePtrVec exp_env_list = env_set_ptr->getExpEnvList();
  std::vector<double> aggregate_inte(exp_env_list[0]->getPeakNum(), 0.0);
  int num_spec = exp_env_list.size();
  int num_peaks = aggregate_inte.size();
  for (int sp_idx = 0; sp_idx < num_spec; sp_idx++) {
    for (int peak_idx = 0; peak_idx < num_peaks; peak_idx++) {
      MatrixPeakPtr peak_ptr = exp_env_list[sp_idx]->getPeakPtr(peak_idx);
      if (peak_ptr != nullptr) {
        aggregate_inte[peak_idx] = aggregate_inte[peak_idx] + peak_ptr->getIntensity();
      }
    }
  }
  return aggregate_inte;
}

std::vector<double> getAggregateEnvelopeMz(EnvSetPtr env_set_ptr) {
  ExpEnvelopePtrVec exp_env_list = env_set_ptr->getExpEnvList();
  std::vector<double> aggregate_mz(exp_env_list[0]->getPeakNum(), 0.0);
  int num_spec = exp_env_list.size();
  int num_peaks = aggregate_mz.size();
  for (int peak_idx = 0; peak_idx < num_peaks; peak_idx++) {
    int counter = 0;
    for (int sp_idx = 0; sp_idx < num_spec; sp_idx++) {
      MatrixPeakPtr peak_ptr = exp_env_list[sp_idx]->getPeakPtr(peak_idx);
      if (peak_ptr != nullptr) {
        aggregate_mz[peak_idx] = aggregate_mz[peak_idx] + peak_ptr->getPosition();
        counter = counter + 1;
      }
    }
    if (counter > 0) {
      aggregate_mz[peak_idx] = aggregate_mz[peak_idx] / counter;
    }
  }
  return aggregate_mz;
}

double calcInteRatio(std::vector<double> &theo_envelope_inte, std::vector<double> &exp_envelope_inte) {
  double theo_sum = 0;
  double obs_sum = 0;
  int refer_idx =
    std::max_element(theo_envelope_inte.begin(), theo_envelope_inte.end()) - theo_envelope_inte.begin();
  theo_sum = theo_sum + theo_envelope_inte[refer_idx];
  obs_sum = obs_sum + exp_envelope_inte[refer_idx];
  if (refer_idx - 1 >= 0) {
    theo_sum = theo_sum + theo_envelope_inte[refer_idx - 1];
    obs_sum = obs_sum + exp_envelope_inte[refer_idx - 1];
  }
  if (refer_idx + 1 < static_cast<int>(theo_envelope_inte.size())) {
    theo_sum = theo_sum + theo_envelope_inte[refer_idx + 1];
    obs_sum = obs_sum + exp_envelope_inte[refer_idx + 1];
  }
  if (theo_sum == 0) {
    return 1.0;
  }
  else {
    return obs_sum / theo_sum;
  }
}

}

}
