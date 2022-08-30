//
// Created by abbash on 8/24/22.
//

#ifndef TOPPIC_ENV_SET_UTIL_HPP
#define TOPPIC_ENV_SET_UTIL_HPP

#include "exp_envelope.hpp"
#include "env_set.hpp"
#include "topfd/feature_detect/spectrum/peak_row.hpp"
#include "topfd/feature_detect/envelope/simple_peak.hpp"
#include "topfd/feature_detect/envelope/seed_envelope.hpp"
#include "topfd/feature_detect/spectrum/peak_matrix.hpp"

namespace toppic {
namespace env_set_util {
  bool check_valid_env_set(PeakMatrix peak_matrix, SeedEnvelope seed_env, double mass_tol, double match_peak_num_tole);
  ExpPeak pick_exp_peak(PeakRow row, SimplePeak seed_peak, double mass_tol);
  ExpEnvelope get_match_exp_env(PeakRow peak_row, SeedEnvelope seed_env, double mass_tol);
  void comp_peak_start_end_idx(PeakMatrix peak_matrix, std::vector<SimplePeak> peak_list, double error_tole);
  void remove_non_match_envs(std::vector<ExpEnvelope> env_list);
  void print_env(ExpEnvelope exp_env, double ratio);
  SeedEnvelope preprocess_env(PeakMatrix peak_matrix, SeedEnvelope seed_env, double mass_tol, bool* valid);
  EnvSet get_env_set_by_three_peaks(PeakMatrix peak_matrix, SeedEnvelope seed_env, double mass_tol, double max_miss_env);
  EnvSet find_env_set(PeakMatrix peak_matrix, SeedEnvelope seed_env, double mass_tol, int start_spec_id, int end_spec_id);
}
}


#endif //TOPPIC_ENV_SET_UTIL_HPP
