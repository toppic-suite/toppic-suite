//
// Created by abbash on 8/29/22.
//

#ifndef TOPPIC_ENV_UTIL_HPP
#define TOPPIC_ENV_UTIL_HPP

#include <vector>
#include "topfd/feature_detect/env_set/env_set.hpp"

namespace toppic {
namespace env_utils {
  double get_mz(double mass, int charge);
  double get_mass(double mz, int charge);
  std::vector<double> getExtendMasses(double mass);

  std::vector<double> get_aggregate_envelopes_mz(EnvSet& env_set);
  std::vector<double> get_aggregate_envelopes_inte(EnvSet& env_set);
  double calcInteRatio_scan(std::vector<double> &theo_envelope_inte, std::vector<double> &exp_envelope_inte);
}
}
#endif //TOPPIC_ENV_UTIL_HPP
