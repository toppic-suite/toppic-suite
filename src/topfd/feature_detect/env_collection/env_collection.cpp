//
// Created by abbash on 7/29/22.
//

#include "env_collection.hpp"
#include <iostream>
#include <valarray>

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
  for (auto env_set : env_set_list_) {
    double cur_weight = 0, cur_weight_mz_error = 0;
    env_set.get_weight_mz_error(&cur_weight, &cur_weight_mz_error);
    weight = weight + cur_weight;
    weight_mz_error = weight_mz_error + cur_weight_mz_error;
  }
  if (weight > 0) {
    double mz_error = weight_mz_error / weight;
    std::cout << "mz error " << mz_error << std::endl;
    seed_env_.shift(mz_error * seed_env_.getCharge());
  }
  else{
   std::cout << "ERROR 0 weight in refine_mono_mass" << std::endl;
  }
}

double toppic::EnvCollection::get_intensity(){
  double inte = 0;
  for (auto env_set : env_set_list_) {
    if (!env_set.isEmpty())
      inte = inte + env_set.comp_intensity();
  }
  return inte;
}

double toppic::EnvCollection::get_min_elution_time(spec_list spectra_list){
  return spectra_list[start_spec_id_].getRt();
}

double toppic::EnvCollection::get_max_elution_time(spec_list spectra_list){
  return spectra_list[end_spec_id_].getRt();
}

double toppic::EnvCollection::get_apex_elution_time(spec_list spectra_list) {
  int ref_spec_id = seed_env_.getSpecId();
  return spectra_list[ref_spec_id].getRt();
}

double toppic::EnvCollection::get_elution_length(spec_list spectra_list) {
  return get_max_elution_time(spectra_list) - get_min_elution_time(spectra_list);
}

double toppic::EnvCollection::remove_matrix_peaks(PeakMatrix peak_matrix) {
  for (auto env_set: env_set_list_) {
    if (!env_set.isEmpty())
      env_set.remove_matrix_peaks(peak_matrix);
  }
}

//void comp_pair_mz_error(){
//  return;
//}
//
//void check_env_len(){
//  return;
//}

