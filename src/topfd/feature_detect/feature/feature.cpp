//
// Created by abbash on 8/30/22.
//

#include "feature.hpp"

toppic::Feature::Feature(EnvCollection &env_coll, PeakMatrix &peak_matrix, fdeep::model &model, int feature_id, double snr){
  SeedEnvelope seed_env = env_coll.getSeedEnv();
  spec_list spectra_list = peak_matrix.get_spectra_list();
  std::vector<std::vector<double>> theo_map = env_coll.get_seed_theo_map(peak_matrix, snr);
  std::vector<double> spectrum_noise_levels = peak_matrix.get_spec_noise_inte();
  double noiseIntensityLevel = std::accumulate(spectrum_noise_levels.begin() + env_coll.getStartSpecId(), spectrum_noise_levels.begin() + env_coll.getEndSpecId(), 0.0);
  int base_spec = env_coll.getBaseSpecID();
  int start_spec = env_coll.getStartSpecId();
  EnvSet env_set = env_coll.get_seed_env_set();

  feature_id_ = feature_id;
  min_scan_ = env_coll.getStartSpecId();
  max_scan_ = env_coll.getEndSpecId();
  min_charge_ = env_coll.getMinCharge();
  max_charge_ = env_coll.getMaxCharge();
  mono_mass_ = seed_env.getMass();
  rep_charge_ = seed_env.getCharge();
  rep_mz_ = seed_env.getPos();
  abundance_ = env_coll.get_intensity(snr, peak_matrix.get_min_inte());
  min_elution_time_ = env_coll.get_min_elution_time(spectra_list)/60.0;
  max_elution_time_ = env_coll.get_max_elution_time(spectra_list)/60.0;
  apex_elution_time_ = env_coll.get_apex_elution_time(spectra_list)/60.0;
  elution_length_ = env_coll.get_elution_length(spectra_list)/60.0;

  percent_matched_peaks_ = component_score::get_matched_peaks_percent(env_set, theo_map);
  intensity_correlation_ = component_score::get_agg_env_corr(env_set);
  top3_correlation_ = component_score::get_3_scan_corr(env_set, base_spec, start_spec);
  even_odd_peak_ratios_ = component_score::get_agg_odd_even_peak_ratio(env_set);
  percent_consec_peaks_ = component_score::get_consecutive_peaks_percent(env_set);
  num_theo_peaks_ = component_score::get_num_theo_peaks(theo_map);
  mz_error_sum_ = component_score::get_mz_errors(env_set);
  envcnn_score_ = env_cnn_score::get_envcnn_score(model, peak_matrix, env_coll, noiseIntensityLevel);
  label_ = 0;
}
