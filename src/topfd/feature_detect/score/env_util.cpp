//
// Created by abbash on 8/29/22.
//

#include "env_util.hpp"
#include <iostream>

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

  std::vector<double> getExtendMasses(double mass) {
    double IM = 1.00235;
    std::vector<double> extend_offsets_;
    for (int i = -1; i <= 1; i++)
      extend_offsets_.push_back(i*IM);
    std::vector<double> result;
    for (double extend_offset : extend_offsets_) {
      double new_mass = mass + extend_offset;
      result.push_back(new_mass);
    }
    return result;
  }

  std::vector<double> get_aggregate_envelopes_inte(EnvSet& env_set){
    std::vector<ExpEnvelope> exp_env_list = env_set.getExpEnvList();
    std::vector<double> aggregate_inte (exp_env_list[0].get_peak_num(), 0.0);
    for (size_t spId = 0; spId < exp_env_list.size(); spId++) {
      for (size_t peakIdx = 0; peakIdx < aggregate_inte.size(); peakIdx++) {
        ExpPeak peak = exp_env_list[spId].get_peak(peakIdx);
        if (!peak.isEmpty())
          aggregate_inte[peakIdx] = aggregate_inte[peakIdx] + peak.getInte();
      }
    }
    return aggregate_inte;
  }

  std::vector<double> get_aggregate_envelopes_mz(EnvSet& env_set){
    std::vector<ExpEnvelope> exp_env_list = env_set.getExpEnvList();
    std::vector<double> aggregate_mz (exp_env_list[0].get_peak_num(), 0.0);
    for (size_t peakIdx = 0; peakIdx < aggregate_mz.size(); peakIdx++) {
      int counter = 0;
      for (size_t spId = 0; spId < exp_env_list.size(); spId++) {
        ExpPeak peak = exp_env_list[spId].get_peak(peakIdx);
        if (!peak.isEmpty()) {
          aggregate_mz[peakIdx] = aggregate_mz[peakIdx] + peak.getPos();
          counter = counter + 1;
        }
      }
      if (counter > 0)
        aggregate_mz[peakIdx] = aggregate_mz[peakIdx]/counter;
    }
    return aggregate_mz;
  }

  double calcInteRatio_scan(std::vector<double> &theo_envelope_inte, std::vector<double> &exp_envelope_inte) {
    double theo_sum = 0;
    double obs_sum = 0;
    int refer_idx = std::max_element(theo_envelope_inte.begin(), theo_envelope_inte.end()) - theo_envelope_inte.begin();
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
