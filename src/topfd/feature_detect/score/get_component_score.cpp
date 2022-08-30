//
// Created by abbash on 8/26/22.
//

#include <valarray>
#include "get_component_score.hpp"

namespace toppic {
namespace component_score {
  std::vector<double> get_theo_envelope_peak_intens(EnvSet env_set){
    std::vector<SimplePeak> peaks = env_set.get_peak_list();
    std::vector<double> theo_peak_intes;
    for (const auto& i: peaks)
      theo_peak_intes.push_back(i.getInte());
    return theo_peak_intes;
  }

  EnvSet get_seed_env_set(EnvCollection env_coll) {
    SeedEnvelope seed_env = env_coll.getSeedEnv();
    std::vector<EnvSet> env_set_list = env_coll.getEnvSetList();
    /// Get Envset of seed envelope
    EnvSet env_set = EnvSet();
    for (const auto& es: env_set_list){
      SeedEnvelope es_seed_env = es.getSeedEnv();
      if (es_seed_env.getCharge() == seed_env.getCharge())
        env_set = EnvSet(es);
    }
    return env_set;
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

  double get_agg_odd_even_peak_ratio(EnvCollection env_coll) {
    EnvSet env_set = get_seed_env_set(env_coll);
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

  double get_agg_env_corr(EnvCollection env_coll) {
    EnvSet env_set = get_seed_env_set(env_coll);
    std::vector<double> theo_inte = get_theo_envelope_peak_intens(env_set);
    std::vector<double> aggregate_inte = env_utils::get_aggregate_envelopes_inte(env_set);
    double max_aggregate_inte = *std::max_element(aggregate_inte.begin(), aggregate_inte.end());
    std::vector<double> normalized_aggregate_inte;
    for (auto inte : aggregate_inte)
      normalized_aggregate_inte.push_back(inte/max_aggregate_inte);
    double corr = utility_functions::pearsonr(theo_inte, normalized_aggregate_inte);
    return corr;
  }

  double get_consecutive_peaks_percent(EnvCollection env_coll) {
    EnvSet env_set = get_seed_env_set(env_coll);
    double total_peaks = 0, positive_peaks = 0;
    std::vector<ExpEnvelope> exp_envs = env_set.getExpEnvList();
    for (auto exp_env : exp_envs) {
      std::vector<ExpPeak> peaks = exp_env.getExpEnvList();
      for (auto p: peaks)
        if (!p.isEmpty())
          total_peaks = total_peaks + 1;
      positive_peaks = positive_peaks + count_max_consecutive_peak_num(exp_env);
    }
    double percent_matched_peaks = positive_peaks/total_peaks;
    return percent_matched_peaks;
  }

  double get_matched_peaks_percent(EnvCollection env_coll) {
    EnvSet env_set = get_seed_env_set(env_coll);
    double total_peaks = 0, positive_peaks = 0;
    std::vector<ExpEnvelope> exp_envs = env_set.getExpEnvList();
    for (auto exp_env : exp_envs) {
      std::vector<ExpPeak> peaks = exp_env.getExpEnvList();
      int sum_val = 0;
      for (auto p: peaks)
        if (!p.isEmpty()) sum_val = sum_val + 1;
      if (sum_val <= 2) continue;
      total_peaks = total_peaks + peaks.size();
      positive_peaks = positive_peaks + sum_val;
    }
    double percent_matched_peaks = -1;
    if (total_peaks > 0) percent_matched_peaks = positive_peaks/total_peaks;
    return percent_matched_peaks;
  }

  double get_3_scan_corr(EnvCollection env_coll) {
    double scan_3_corr = 0;
    EnvSet env_set = get_seed_env_set(env_coll);
    std::vector<ExpEnvelope> exp_envs = env_set.getExpEnvList();
    std::vector<ExpEnvelope> shortlisted_spec;
    for (const auto& exp_env : exp_envs){
      SeedEnvelope seed_env = env_set.getSeedEnv();
      if (seed_env.getSpecId() - 1 <= exp_env.getSpecId() <=  seed_env.getSpecId() + 1)
        shortlisted_spec.push_back(exp_env);
    }
    std::vector<double> exp_peaks_base_spec_1;
    std::vector<double> exp_peaks_base_spec_2;
    std::vector<double> exp_peaks_base_spec_3;
    if (shortlisted_spec.size() == 2) {
      std::vector<ExpPeak> peaks = shortlisted_spec[0].getExpEnvList();
      for (auto peak : peaks)
        if(!peak.isEmpty()) exp_peaks_base_spec_1.push_back(peak.getInte());
      peaks = shortlisted_spec[1].getExpEnvList();
      for (auto peak : peaks)
        if(!peak.isEmpty()) exp_peaks_base_spec_2.push_back(peak.getInte());
      scan_3_corr = utility_functions::pearsonr(exp_peaks_base_spec_1, exp_peaks_base_spec_2);
    }
    if (shortlisted_spec.size() == 3) {
      std::vector<ExpPeak> peaks = shortlisted_spec[0].getExpEnvList();
      for (auto peak : peaks)
        if(!peak.isEmpty()) exp_peaks_base_spec_1.push_back(peak.getInte());
      peaks = shortlisted_spec[1].getExpEnvList();
      for (auto peak : peaks)
        if(!peak.isEmpty()) exp_peaks_base_spec_2.push_back(peak.getInte());
      peaks = shortlisted_spec[2].getExpEnvList();
      for (auto peak : peaks)
        if(!peak.isEmpty()) exp_peaks_base_spec_3.push_back(peak.getInte());
      double corr_sp_12 = utility_functions::pearsonr(exp_peaks_base_spec_1, exp_peaks_base_spec_2);
      double corr_sp_13 = utility_functions::pearsonr(exp_peaks_base_spec_1, exp_peaks_base_spec_3);
      double corr_sp_23 = utility_functions::pearsonr(exp_peaks_base_spec_2, exp_peaks_base_spec_3);
      scan_3_corr = (corr_sp_12 + corr_sp_13 + corr_sp_23)/3.0;
    }
    return scan_3_corr;
  }

  double get_rt_range(EnvCollection env_coll){
    return (env_coll.getEndSpecId() - env_coll.getStartSpecId());
//    return (env_coll. get_max_elution_time() - env_coll.get_min_elution_time());
  }

  double get_charge_range(EnvCollection env_coll) {
    return (env_coll.getMaxCharge() - env_coll.getMinCharge());
  }

}
}
