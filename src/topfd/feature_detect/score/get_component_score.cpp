//
// Created by abbash on 8/26/22.
//

#include <valarray>
#include "get_component_score.hpp"

namespace toppic {
namespace component_score {
  std::vector<double> get_theo_envelope_peak_intens(EnvSet& env_set){
    std::vector<SimplePeak> peaks = env_set.get_peak_list();
    std::vector<double> theo_peak_intes;
    for (const auto& i: peaks)
      theo_peak_intes.push_back(i.getInte());
    return theo_peak_intes;
  }

  double count_max_consecutive_peak_num(ExpEnvelope exp_env) {
    int n = 0, max_consecutive_peak_num = 0;
    std::vector<ExpPeak> peaks = exp_env.getExpEnvList();
    for (auto peak : peaks){
      if (!peak.isEmpty()){
        n = n + 1;
        if (n > max_consecutive_peak_num)
          max_consecutive_peak_num = 0;
      }
      else
        n = 0;
    }
    return max_consecutive_peak_num;
  }

  double get_agg_odd_even_peak_ratio(EnvSet& env_set) {
    std::vector<double> theo_inte = get_theo_envelope_peak_intens(env_set);
    std::vector<double> aggregate_inte = env_utils::get_aggregate_envelopes_inte(env_set);
    double max_aggregate_inte = *std::max_element(aggregate_inte.begin(), aggregate_inte.end());
    std::vector<double> normalized_aggregate_inte;
    for (auto inte : aggregate_inte)
      normalized_aggregate_inte.push_back(inte/max_aggregate_inte);
    double sum_even_peaks = 0, sum_even_peaks_theo = 0, sum_odd_peaks = 0, sum_odd_peaks_theo = 0;
    for (int peak_idx = 0; peak_idx < aggregate_inte.size(); peak_idx++) {
      if (peak_idx%2 == 0) {
        sum_even_peaks = sum_even_peaks + normalized_aggregate_inte[peak_idx];
        sum_even_peaks_theo = sum_even_peaks_theo + theo_inte[peak_idx];
      }
      else {
        sum_odd_peaks = sum_odd_peaks + normalized_aggregate_inte[peak_idx];
        sum_odd_peaks_theo = sum_odd_peaks_theo + theo_inte[peak_idx];
      }
    }
    double inte_ratio = (sum_even_peaks/sum_even_peaks_theo) / (sum_odd_peaks/sum_odd_peaks_theo);
    if (inte_ratio <= 0)
      inte_ratio = 1;
    return std::log10(inte_ratio);
  }

  double get_agg_env_corr(EnvSet& env_set) {
    std::vector<double> theo_inte = get_theo_envelope_peak_intens(env_set);
    std::vector<double> aggregate_inte = env_utils::get_aggregate_envelopes_inte(env_set);
    double max_aggregate_inte = *std::max_element(aggregate_inte.begin(), aggregate_inte.end());
    std::vector<double> normalized_aggregate_inte;
    for (auto inte : aggregate_inte)
      normalized_aggregate_inte.push_back(inte/max_aggregate_inte);
    double corr = utility_functions::pearsonr(theo_inte, normalized_aggregate_inte);
    return corr;
  }

  double get_mz_errors(EnvSet& env_set) {
    std::vector<double> theo_dis = env_set.get_theo_distribution_mz();
    std::vector<ExpEnvelope> exp_envs = env_set.getExpEnvList();
    double error_sum = 0;
    for (auto & exp_env : exp_envs) {
      std::vector<ExpPeak> peaks = exp_env.getExpEnvList();
      for (int peak_idx = 0; peak_idx < peaks.size(); peak_idx++) {
        ExpPeak peak = peaks[peak_idx];
        if (!peak.isEmpty()) {
          double cur_err = std::abs(peak.getPos() - theo_dis[peak_idx]);
          error_sum = error_sum + cur_err;
        }
      }
    }
//    std::cout << "Error SUM: " << error_sum << std::endl;
    return error_sum;
  }

  double count_max_consecutive_peak_num(std::vector<ExpPeak>& peaks) {
    int n = 0;
    int max_consecutive_peak_num = 0;
    for (auto peak : peaks){
      if (!peak.isEmpty()) {
        n = n + 1;
        if (n > max_consecutive_peak_num)
          max_consecutive_peak_num = n;
      }
      else
        n = 0;
    }
    return max_consecutive_peak_num;
  }

  double get_consecutive_peaks_percent(EnvSet& env_set) {
    double total_peaks = 0, positive_peaks = 0;
    std::vector<ExpEnvelope> exp_envs = env_set.getExpEnvList();
    for (auto & exp_env : exp_envs) {
      std::vector<ExpPeak> peaks = exp_env.getExpEnvList();
      for (auto p : peaks)
        if (!p.isEmpty())
          total_peaks = total_peaks + 1;
      positive_peaks  = positive_peaks + count_max_consecutive_peak_num(peaks);
    }
    double percent_matched_peaks = positive_peaks/total_peaks;
    return percent_matched_peaks;
  }

  double get_matched_peaks_percent(EnvSet& env_set, std::vector<std::vector<double>> theo_map) {
    double total_peaks = 0, positive_peaks = 0;
    std::vector<ExpEnvelope> exp_envs = env_set.getExpEnvList();
    for (int i = 0; i < exp_envs.size(); i++) {
      std::vector<ExpPeak> peaks = exp_envs[i].getExpEnvList();
      std::vector<double> scalled_theo_env = theo_map[i];
      for (int peak_id = 0; peak_id < scalled_theo_env.size(); peak_id++){
        double peak_inte = scalled_theo_env[peak_id];
        if (peak_inte > 0){
          total_peaks = total_peaks + 1;
          if (!peaks[peak_id].isEmpty())
            positive_peaks = positive_peaks + 1;
        }
      }
//      std::cout << i << ", " << exp_envs[i].getSpecId() << ", " << total_peaks << ", " << positive_peaks << ", " << positive_peaks/total_peaks << std::endl;
    }
    double percent_matched_peaks = -1;
    if (total_peaks > 0) percent_matched_peaks = positive_peaks/total_peaks;
    return percent_matched_peaks;
  }

  double get_num_theo_peaks(std::vector<std::vector<double>>& theo_map) {
    double total_peaks = 0;
    for (auto & scan_intes : theo_map)
      for (double peak_inte : scan_intes)
        if (peak_inte > 0)
          total_peaks = total_peaks + 1;
    return total_peaks;
  }

  double get_3_scan_corr(EnvSet& env_set, int base_spec, int start_spec) {
    double scan_3_corr;
    std::vector<ExpEnvelope> exp_envs = env_set.getExpEnvList();
    base_spec = std::max(base_spec - start_spec, 0);

//    std::cout << base_spec << ", " << exp_envs.size() << ", " << env_coll.getBaseSpecID() << ", " <<  env_coll.getStartSpecId() << std::endl;
    std::vector<double> data_sp = exp_envs[base_spec].get_inte_list();
    std::vector<double> data_sp_minus_1 (data_sp.size(), 0.0);
    std::vector<double> data_sp_plus_1  (data_sp.size(), 0.0);

    if (base_spec - 1 > 0)
      data_sp_minus_1 = exp_envs[base_spec-1].get_inte_list();
    if (base_spec + 1 < exp_envs.size())
      data_sp_plus_1 = exp_envs[base_spec+1].get_inte_list();

//    for (int i = 0; i < data_sp.size(); i++)
//      std::cout << i << ", " << data_sp_minus_1[i] << ", " << data_sp[i] << ", " << data_sp_plus_1[i] << std::endl;

    double sp_sum = std::accumulate(data_sp.begin(), data_sp.end(), 0.0);
    double sp_minus_1_sum = std::accumulate(data_sp_minus_1.begin(), data_sp_minus_1.end(), 0.0);
    double sp_plus_1_sum = std::accumulate(data_sp_plus_1.begin(), data_sp_plus_1.end(), 0.0);
//    std::cout << "SUM: " << sp_minus_1_sum << ", " << sp_sum << ", " << sp_plus_1_sum << std::endl;

    if (sp_sum > 0 and sp_minus_1_sum > 0 and sp_plus_1_sum > 0){
      double corr_sp_12 = utility_functions::pearsonr(data_sp, data_sp_minus_1);
      double corr_sp_13 = utility_functions::pearsonr(data_sp, data_sp_plus_1);
      double corr_sp_23 = utility_functions::pearsonr(data_sp_minus_1, data_sp_plus_1);
      scan_3_corr = (corr_sp_12 + corr_sp_13 + corr_sp_23)/3.0;
    }
    else if (sp_sum > 0 and sp_minus_1_sum > 0 and sp_plus_1_sum == 0){
      scan_3_corr = utility_functions::pearsonr(data_sp, data_sp_minus_1);
    }
    else if (sp_sum > 0 and sp_minus_1_sum == 0 and sp_plus_1_sum > 0){
      scan_3_corr = utility_functions::pearsonr(data_sp, data_sp_plus_1);
    }
    else if (sp_sum == 0 and sp_minus_1_sum > 0 and sp_plus_1_sum > 0){
      scan_3_corr = utility_functions::pearsonr(data_sp_minus_1, data_sp_plus_1);
    }
    else
      scan_3_corr = 0;
//    std::cout << "Top 3 scans correlation: " << scan_3_corr << std::endl;
    return scan_3_corr;
  }

  //  double get_3_scan_corr(EnvCollection& env_coll) {
//    double scan_3_corr = 0;
//    EnvSet env_set = get_seed_env_set(env_coll);
//    std::vector<ExpEnvelope> exp_envs = env_set.getExpEnvList();
//
//    std::vector<ExpEnvelope> shortlisted_spec;
//    for (const auto& exp_env : exp_envs){
//      SeedEnvelope seed_env = env_set.getSeedEnv();
//      if (seed_env.getSpecId() - 1 <= exp_env.getSpecId() <=  seed_env.getSpecId() + 1)
//        shortlisted_spec.push_back(exp_env);
//    }
//    std::vector<double> exp_peaks_base_spec_1;
//    std::vector<double> exp_peaks_base_spec_2;
//    std::vector<double> exp_peaks_base_spec_3;
//    if (shortlisted_spec.size() == 2) {
//      std::vector<ExpPeak> peaks = shortlisted_spec[0].getExpEnvList();
//      for (auto peak : peaks)
//        if(!peak.isEmpty()) exp_peaks_base_spec_1.push_back(peak.getInte());
//      peaks = shortlisted_spec[1].getExpEnvList();
//      for (auto peak : peaks)
//        if(!peak.isEmpty()) exp_peaks_base_spec_2.push_back(peak.getInte());
//      scan_3_corr = utility_functions::pearsonr(exp_peaks_base_spec_1, exp_peaks_base_spec_2);
//    }
//    if (shortlisted_spec.size() == 3) {
//      std::vector<ExpPeak> peaks = shortlisted_spec[0].getExpEnvList();
//      for (auto peak : peaks)
//        if(!peak.isEmpty()) exp_peaks_base_spec_1.push_back(peak.getInte());
//      peaks = shortlisted_spec[1].getExpEnvList();
//      for (auto peak : peaks)
//        if(!peak.isEmpty()) exp_peaks_base_spec_2.push_back(peak.getInte());
//      peaks = shortlisted_spec[2].getExpEnvList();
//      for (auto peak : peaks)
//        if(!peak.isEmpty()) exp_peaks_base_spec_3.push_back(peak.getInte());
//      double corr_sp_12 = utility_functions::pearsonr(exp_peaks_base_spec_1, exp_peaks_base_spec_2);
//      double corr_sp_13 = utility_functions::pearsonr(exp_peaks_base_spec_1, exp_peaks_base_spec_3);
//      double corr_sp_23 = utility_functions::pearsonr(exp_peaks_base_spec_2, exp_peaks_base_spec_3);
//      scan_3_corr = (corr_sp_12 + corr_sp_13 + corr_sp_23)/3.0;
//    }
//    return scan_3_corr;
//  }

  double get_rt_range(EnvCollection& env_coll){
    return (env_coll.getEndSpecId() - env_coll.getStartSpecId());
//    return (env_coll. get_max_elution_time() - env_coll.get_min_elution_time());
  }

  double get_charge_range(EnvCollection env_coll) {
    return (env_coll.getMaxCharge() - env_coll.getMinCharge());
  }

}
}
