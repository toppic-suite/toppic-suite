//
// Created by abbash on 9/5/22.
//

#include "evaluate_envelope.hpp"
#include <iostream>

bool toppic::evaluate_envelope::preprocess_env(PeakMatrix &peak_matrix, SeedEnvelope &env, double mass_tol, double corr_tol, bool valid) {
  double min_mz = peak_matrix.get_min_mz() - mass_tol;
  double max_mz = peak_matrix.get_max_mz() + mass_tol;
  env.rm_peaks(min_mz, max_mz);
  env_set_util::comp_peak_start_end_idx(peak_matrix, env, mass_tol);
  valid = evaluate_envelope(peak_matrix, env, mass_tol, corr_tol);
  if (env.getSpecId() >= peak_matrix.get_spec_num()) {
    std::cout << "spec id " << env.getSpecId() << " is out of range!" << std::endl;
    valid = false;
  }
  return valid;
}

bool toppic::evaluate_envelope::evaluate_envelope(PeakMatrix &peak_matrix, SeedEnvelope &seed_envelope, double mass_tol, double corr_tol){
  double snr = 3.0;
  double noise_inte = peak_matrix.get_min_inte();
//  PeakRow row = peak_matrix.get_row(seed_envelope.getSpecId());
//  ExpEnvelope exp_env = env_set_util::get_match_exp_env(row, seed_envelope, mass_tol);
  ExpEnvelope exp_env = env_set_util::get_match_exp_env(peak_matrix, seed_envelope, seed_envelope.getSpecId(), mass_tol);
  std::vector<double> experimental_envelope_mass = exp_env.get_pos_list();
  std::vector<double> seed_envelope_mass = seed_envelope.get_pos_list();
  std::vector<double> experimental_envelope_inte = exp_env.get_inte_list();
  std::vector<double> seed_envelope_inte = seed_envelope.get_inte_list();
  int num_env_peaks = seed_envelope_inte.size();
  double inte_ratio = env_utils::calcInteRatio_scan(seed_envelope_inte, experimental_envelope_inte);
  std::vector<double> scaled_theo_inte;

  int i;
  for (i = 0; i < num_env_peaks; i++) {
    double scaled_inte = inte_ratio * seed_envelope_inte[i];
    if (scaled_inte < snr * noise_inte)
      break;
    scaled_theo_inte.push_back(scaled_inte);
  }
  std::vector<SimplePeak> seed_envelope_peaks = seed_envelope.getPeakList();
  for (int j = num_env_peaks-1; j >= i; j--) {
    seed_envelope_peaks.erase(seed_envelope_peaks.begin() + j);
    experimental_envelope_inte.erase(experimental_envelope_inte.begin() + j);
  }
  seed_envelope.setPeakList(seed_envelope_peaks);
  if (!test_charge_state(seed_envelope.getCharge(), scaled_theo_inte)) return false;
  return evaluate_envelope_pair(experimental_envelope_inte, scaled_theo_inte, corr_tol);
}

bool toppic::evaluate_envelope::test_charge_state(int charge, std::vector<double> &seed_envelope_inte) {
  if ((charge == 1 || charge == 2) && seed_envelope_inte.size() < 2) return false;
  if (charge > 2 && charge < 15 && seed_envelope_inte.size() < 3) return false;
  if (charge >= 15 && seed_envelope_inte.size() < 5) return false;
  return true;
}

bool toppic::evaluate_envelope::evaluate_envelope_pair(std::vector<double> &experimental_envelope_inte, std::vector<double> &theo_inte, double tol) {
  int num_theo_peak = 0;
  for (auto i : theo_inte) if (i > 0) num_theo_peak++;
  if (num_theo_peak == 0) return false;
  double corr = toppic::utility_functions::pearsonr(experimental_envelope_inte, theo_inte);
//  std::cout << "Correlation of Seed Envelope: " << corr << " and number of scalled theo peaks " << num_theo_peak << std::endl;
  if (corr < tol) return false;
  return true;
}


