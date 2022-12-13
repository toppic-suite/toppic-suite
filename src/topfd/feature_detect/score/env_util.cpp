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

#include "env_util.hpp"

namespace toppic {
  namespace env_utils {
    double get_mz(double mass, int charge) {
      double proton = 1.00727;
      return (mass + (charge * proton)) / charge;
    }

    double get_mass(double mz, int charge) {
      double proton = 1.00727;
      return (mz * charge) - (charge * proton);
    }

    std::vector<double> get_aggregate_envelopes_inte(EnvSet &env_set) {
      std::vector<ExpEnvelope> exp_env_list = env_set.getExpEnvList();
      std::vector<double> aggregate_inte(exp_env_list[0].get_peak_num(), 0.0);
      int num_spec = exp_env_list.size();
      int num_peaks = aggregate_inte.size();
      for (int spId = 0; spId < num_spec; spId++) {
        for (int peakIdx = 0; peakIdx < num_peaks; peakIdx++) {
          ExpPeak peak = exp_env_list[spId].get_peak(peakIdx);
          if (!peak.isEmpty())
            aggregate_inte[peakIdx] = aggregate_inte[peakIdx] + peak.getInte();
        }
      }
      return aggregate_inte;
    }

    std::vector<double> get_aggregate_envelopes_mz(EnvSet &env_set) {
      std::vector<ExpEnvelope> exp_env_list = env_set.getExpEnvList();
      std::vector<double> aggregate_mz(exp_env_list[0].get_peak_num(), 0.0);
      int num_spec = exp_env_list.size();
      int num_peaks = aggregate_mz.size();
      for (int peakIdx = 0; peakIdx < num_peaks; peakIdx++) {
        int counter = 0;
        for (int spId = 0; spId < num_spec; spId++) {
          ExpPeak peak = exp_env_list[spId].get_peak(peakIdx);
          if (!peak.isEmpty()) {
            aggregate_mz[peakIdx] = aggregate_mz[peakIdx] + peak.getPos();
            counter = counter + 1;
          }
        }
        if (counter > 0)
          aggregate_mz[peakIdx] = aggregate_mz[peakIdx] / counter;
      }
      return aggregate_mz;
    }

    double calcInteRatio_scan(std::vector<double> &theo_envelope_inte, std::vector<double> &exp_envelope_inte) {
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
      if (theo_sum == 0)
        return 1.0;
      else
        return obs_sum / theo_sum;
    }
  }
}
