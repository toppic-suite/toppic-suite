//
// Created by abbash on 7/29/22.
//

#include "env_collection.hpp"
#include "topfd/feature_detect/score/env_util.hpp"
#include <iostream>
#include <valarray>

toppic::EnvCollection::EnvCollection() {
  seed_env_ = SeedEnvelope();
  min_charge_ = -1;
  max_charge_ = -1;
  start_spec_id_ = -1;
  end_spec_id_ = -1;
}

toppic::EnvCollection::EnvCollection(const EnvCollection &ec) {
  seed_env_ = ec.seed_env_;
  min_charge_ = ec.min_charge_;
  max_charge_ = ec.max_charge_;
  start_spec_id_ = ec.start_spec_id_;
  end_spec_id_ = ec.end_spec_id_;
  for (auto & i : ec.env_set_list_)
    env_set_list_.push_back(i);
  for (auto & i : ec.exp_inte_sum_list_)
    exp_inte_sum_list_.push_back(i);
}

toppic::EnvCollection::EnvCollection(const SeedEnvelope &env, const std::vector<EnvSet> &env_set_list, int min_charge, int max_charge, int start_spec_id, int end_spec_id){
  seed_env_ = env;
  for (auto & i : env_set_list) env_set_list_.push_back(i);
  min_charge_ = min_charge;
  max_charge_ = max_charge;
  start_spec_id_ = start_spec_id;
  end_spec_id_ = end_spec_id;
  exp_inte_sum_list_ = comp_exp_inte_sum_list();
}

std::vector<double> toppic::EnvCollection::comp_exp_inte_sum_list() {
  int peak_num = seed_env_.get_peak_num();
  std::vector<double> sum_list (peak_num, 0);
  for (auto & env_set : env_set_list_) {
    if (env_set.isEmpty())
      continue;
    std::vector<double> cur_sum_list = env_set.comp_exp_inte_sum_list();
    if (peak_num != cur_sum_list.size())
      std::cout << "peak number " << peak_num << " cur list len " << cur_sum_list.size();
    for (int i = 0; i < peak_num; i++)
      sum_list[i] = sum_list[i] + cur_sum_list[i];
  }
  return sum_list;
}

double toppic::EnvCollection::comp_correlation() {
  std::vector<double> theo_dist = seed_env_.get_inte_list();
  std::vector<double> exp_dist = exp_inte_sum_list_;
  if (theo_dist.size() != exp_dist.size()) {
    std::cout << theo_dist.size() << ", " << exp_dist.size() << std::endl;
    return 0;
  }
  double correlation = utility_functions::pearsonr(theo_dist, exp_dist);
  return correlation;
}

double toppic::EnvCollection::comp_odd_even_log_ratio(){
  std::vector<double> ref_inte_list = seed_env_.get_inte_list();
  int peak_num = ref_inte_list.size();

  double odd_ref_inte = 0;
  double odd_exp_inte = 0;
  for (int i = 0; i <= peak_num; i = i + 2) {
    odd_ref_inte = odd_ref_inte + ref_inte_list[i];
    odd_exp_inte = odd_exp_inte + exp_inte_sum_list_[i];
  }
  double odd_ratio = odd_exp_inte / odd_ref_inte;

  double even_ref_inte = 0, even_exp_inte = 0;
  for (int i = 1; i <= peak_num; i = i + 2) {
    even_ref_inte = even_ref_inte + ref_inte_list[i];
    even_exp_inte = even_exp_inte + exp_inte_sum_list_[i];
  }
  double even_ratio = even_exp_inte/even_ref_inte;

  if (even_ratio == 0 || odd_ratio == 0)
    return 5;
  double log_ratio = std::abs(log10(odd_ratio/even_ratio));
  return log_ratio;
}

void toppic::EnvCollection::refine_mono_mass(){
  double weight = 0;
  double weight_mz_error = 0;
  for (auto &env_set : env_set_list_) {
    double cur_weight = 0, cur_weight_mz_error = 0;
    env_set.get_weight_mz_error(&cur_weight, &cur_weight_mz_error);
    weight = weight + cur_weight;
    weight_mz_error = weight_mz_error + cur_weight_mz_error;
  }
  if (weight > 0) {
    double mz_error = weight_mz_error / weight;
    seed_env_.shift(mz_error * seed_env_.getCharge());
  }
  else{
   std::cout << "ERROR 0 weight in refine_mono_mass" << std::endl;
  }
}

double toppic::EnvCollection::get_intensity(double snr, double noise_inte){
  double inte = 0;
  for (auto env_set : env_set_list_) {
    double tmp_inte = env_set.comp_intensity(snr, noise_inte);
//    std::cout << env_set.getCharge() << ", " << tmp_inte << std::endl;
    inte = inte + tmp_inte;
  }
  return inte;
}

double toppic::EnvCollection::get_min_elution_time(spec_list &spectra_list){
  return spectra_list[start_spec_id_].getRt();
}

double toppic::EnvCollection::get_max_elution_time(spec_list &spectra_list){
  return spectra_list[end_spec_id_].getRt();
}

double toppic::EnvCollection::get_apex_elution_time(spec_list &spectra_list) {
  int ref_spec_id = seed_env_.getSpecId();
  return spectra_list[ref_spec_id].getRt();
}

double toppic::EnvCollection::get_elution_length(spec_list &spectra_list) {
  return get_max_elution_time(spectra_list) - get_min_elution_time(spectra_list);
}

void toppic::EnvCollection::remove_matrix_peaks(PeakMatrix &peak_matrix) {
  for (auto env_set: env_set_list_) {
    if (!env_set.isEmpty())
      env_set.remove_matrix_peaks(peak_matrix);
  }
}

void toppic::EnvCollection::remove_peak_data(PeakMatrix &peak_matrix) {
  for (auto env_set: env_set_list_) {
    if (!env_set.isEmpty())
      env_set.remove_peak_data(peak_matrix);
  }
}

std::vector<std::vector<double>> toppic::EnvCollection::get_seed_theo_map(PeakMatrix &peak_matrix, double snr) {
  double noise_inte = peak_matrix.get_min_inte();
  SeedEnvelope seed_env = this->seed_env_;
  std::vector<EnvSet> env_set_list = this->env_set_list_;
  EnvSet env_set = EnvSet();
  for (auto &es: env_set_list) {
    SeedEnvelope es_seed_env = es.getSeedEnv();
    if (es_seed_env.getCharge() == seed_env.getCharge())
      env_set = EnvSet(es);
  }
  std::vector<std::vector<double>> map = env_set.get_map(snr, noise_inte);
  return map;
}

toppic::EnvSet toppic::EnvCollection::get_seed_env_set() {
  EnvSet env_set = EnvSet();
  for (const auto &es: env_set_list_) {
    SeedEnvelope es_seed_env = es.getSeedEnv();
    if (es_seed_env.getCharge() == seed_env_.getCharge())
      env_set = EnvSet(es);
  }
  return env_set;
}
