//
// Created by abbash on 8/26/22.
//

#ifndef TOPPIC_UTILITY_FUNCTIONS_HPP
#define TOPPIC_UTILITY_FUNCTIONS_HPP

#include <vector>
#include <iostream>
#include "topfd/feature_detect/envelope/seed_envelope.hpp"
#include "topfd/feature_detect/env_set/env_set.hpp"
#include "topfd/feature_detect/env_collection/env_collection.hpp"
#include "topfd/feature_detect/envelope/evaluate_envelope.hpp"


namespace toppic {
namespace utility_functions {
  double pearsonr(std::vector<double> X, std::vector<double> Y);
  std::vector<double> findLocalMinima(std::vector<double> arr);
  std::vector<double> findLocalMaxima(std::vector<double> arr);
  SeedEnvelope test_half_charge_state(PeakMatrix& peak_matrix, SeedEnvelope& env, EnvSet& top_peak_env_set, double even_odd_peak_ratios, double mass_tole);
//  bool check_in_existing_features(PeakMatrix& peakMatrix, EnvCollection& env_coll, std::vector<EnvCollection> env_coll_list, double match_envelope_tolerance, double time_tol);
}
}


#endif //TOPPIC_UTILITY_FUNCTIONS_HPP
