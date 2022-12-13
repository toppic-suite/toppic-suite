//Copyright (c) 2014 - 2022, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#include "get_component_score.hpp"

namespace toppic {
  namespace component_score {

    std::vector<double> get_theo_envelope_peak_intens(EnvSet &env_set) {
      std::vector<SimplePeak> peaks = env_set.get_peak_list();
      std::vector<double> theo_peak_intes;
      for (const auto &i: peaks)
        theo_peak_intes.push_back(i.getInte());
      return theo_peak_intes;
    }

    double get_agg_odd_even_peak_ratio(EnvSet &env_set) {
      std::vector<double> theo_inte = get_theo_envelope_peak_intens(env_set);
      std::vector<double> aggregate_inte = env_utils::get_aggregate_envelopes_inte(env_set);
      double max_aggregate_inte = *std::max_element(aggregate_inte.begin(), aggregate_inte.end());
      int num_peaks = aggregate_inte.size();
      std::vector<double> normalized_aggregate_inte;
      for (auto inte: aggregate_inte)
        normalized_aggregate_inte.push_back(inte / max_aggregate_inte);
      double sum_even_peaks = 0, sum_even_peaks_theo = 0, sum_odd_peaks = 0, sum_odd_peaks_theo = 0;
      for (int peak_idx = 0; peak_idx < num_peaks; peak_idx++) {
        if (peak_idx % 2 == 0) {
          sum_even_peaks = sum_even_peaks + normalized_aggregate_inte[peak_idx];
          sum_even_peaks_theo = sum_even_peaks_theo + theo_inte[peak_idx];
        } else {
          sum_odd_peaks = sum_odd_peaks + normalized_aggregate_inte[peak_idx];
          sum_odd_peaks_theo = sum_odd_peaks_theo + theo_inte[peak_idx];
        }
      }
      double inte_ratio = (sum_even_peaks / sum_even_peaks_theo) / (sum_odd_peaks / sum_odd_peaks_theo);
      if (inte_ratio <= 0)
        inte_ratio = 1;
      return std::log10(inte_ratio);
    }

    double get_agg_env_corr(EnvSet &env_set) {
      std::vector<double> theo_inte = get_theo_envelope_peak_intens(env_set);
      std::vector<double> aggregate_inte = env_utils::get_aggregate_envelopes_inte(env_set);
      double max_aggregate_inte = *std::max_element(aggregate_inte.begin(), aggregate_inte.end());
      std::vector<double> normalized_aggregate_inte;
      for (auto inte: aggregate_inte)
        normalized_aggregate_inte.push_back(inte / max_aggregate_inte);
      double corr = utility_functions::pearsonr(theo_inte, normalized_aggregate_inte);
      return corr;
    }

    double get_mz_errors(EnvSet &env_set) {
      std::vector<double> theo_dis = env_set.get_theo_distribution_mz();
      std::vector<ExpEnvelope> exp_envs = env_set.getExpEnvList();
      double error_sum = 0;
      for (auto &exp_env: exp_envs) {
        std::vector<ExpPeak> peaks = exp_env.getExpEnvList();
        int num_peaks = peaks.size();
        for (int peak_idx = 0; peak_idx < num_peaks; peak_idx++) {
          ExpPeak peak = peaks[peak_idx];
          if (!peak.isEmpty()) {
            double cur_err = std::abs(peak.getPos() - theo_dis[peak_idx]);
            error_sum = error_sum + cur_err;
          }
        }
      }
      return error_sum;
    }

    double count_max_consecutive_peak_num(std::vector<ExpPeak> &peaks) {
      int n = 0;
      int max_consecutive_peak_num = 0;
      for (auto peak: peaks) {
        if (!peak.isEmpty()) {
          n = n + 1;
          if (n > max_consecutive_peak_num)
            max_consecutive_peak_num = n;
        } else
          n = 0;
      }
      return max_consecutive_peak_num;
    }

    double get_consecutive_peaks_percent(EnvSet &env_set) {
      double total_peaks = 0, positive_peaks = 0;
      std::vector<ExpEnvelope> exp_envs = env_set.getExpEnvList();
      for (auto &exp_env: exp_envs) {
        std::vector<ExpPeak> peaks = exp_env.getExpEnvList();
        for (auto p: peaks)
          if (!p.isEmpty())
            total_peaks = total_peaks + 1;
        positive_peaks = positive_peaks + count_max_consecutive_peak_num(peaks);
      }
      double percent_matched_peaks = positive_peaks / total_peaks;
      return percent_matched_peaks;
    }

    double get_matched_peaks_percent(EnvSet &env_set, std::vector<std::vector<double>> theo_map) {
      double total_peaks = 0, positive_peaks = 0;
      std::vector<ExpEnvelope> exp_envs = env_set.getExpEnvList();
      int num_exp_envs = exp_envs.size();
      for (int i = 0; i < num_exp_envs; i++) {
        std::vector<ExpPeak> peaks = exp_envs[i].getExpEnvList();
        std::vector<double> scalled_theo_env = theo_map[i];
        int num_peaks = scalled_theo_env.size();
        for (int peak_id = 0; peak_id < num_peaks; peak_id++) {
          double peak_inte = scalled_theo_env[peak_id];
          if (peak_inte > 0) {
            total_peaks = total_peaks + 1;
            if (!peaks[peak_id].isEmpty())
              positive_peaks = positive_peaks + 1;
          }
        }
      }
      double percent_matched_peaks = -1;
      if (total_peaks > 0) percent_matched_peaks = positive_peaks / total_peaks;
      return percent_matched_peaks;
    }

    int get_num_theo_peaks(std::vector<std::vector<double>> &theo_map) {
      int total_peaks = 0;
      for (auto &scan_intes: theo_map)
        for (double peak_inte: scan_intes)
          if (peak_inte > 0)
            total_peaks = total_peaks + 1;
      return total_peaks;
    }

    double get_3_scan_corr(EnvSet &env_set, int base_spec, int start_spec) {
      double scan_3_corr;
      std::vector<ExpEnvelope> exp_envs = env_set.getExpEnvList();
      base_spec = std::max(base_spec - start_spec, 0);
      std::vector<double> data_sp = exp_envs[base_spec].get_inte_list();
      std::vector<double> data_sp_minus_1(data_sp.size(), 0.0);
      std::vector<double> data_sp_plus_1(data_sp.size(), 0.0);

      if (base_spec - 1 > 0)
        data_sp_minus_1 = exp_envs[base_spec - 1].get_inte_list();
      if (base_spec + 1 < static_cast<int>(exp_envs.size()))
        data_sp_plus_1 = exp_envs[base_spec + 1].get_inte_list();
      double sp_sum = std::accumulate(data_sp.begin(), data_sp.end(), 0.0);
      double sp_minus_1_sum = std::accumulate(data_sp_minus_1.begin(), data_sp_minus_1.end(), 0.0);
      double sp_plus_1_sum = std::accumulate(data_sp_plus_1.begin(), data_sp_plus_1.end(), 0.0);
      if (sp_sum > 0 and sp_minus_1_sum > 0 and sp_plus_1_sum > 0) {
        double corr_sp_12 = utility_functions::pearsonr(data_sp, data_sp_minus_1);
        double corr_sp_13 = utility_functions::pearsonr(data_sp, data_sp_plus_1);
        double corr_sp_23 = utility_functions::pearsonr(data_sp_minus_1, data_sp_plus_1);
        scan_3_corr = (corr_sp_12 + corr_sp_13 + corr_sp_23) / 3.0;
      } else if (sp_sum > 0 and sp_minus_1_sum > 0 and sp_plus_1_sum == 0) {
        scan_3_corr = utility_functions::pearsonr(data_sp, data_sp_minus_1);
      } else if (sp_sum > 0 and sp_minus_1_sum == 0 and sp_plus_1_sum > 0) {
        scan_3_corr = utility_functions::pearsonr(data_sp, data_sp_plus_1);
      } else if (sp_sum == 0 and sp_minus_1_sum > 0 and sp_plus_1_sum > 0) {
        scan_3_corr = utility_functions::pearsonr(data_sp_minus_1, data_sp_plus_1);
      } else
        scan_3_corr = 0;
      return scan_3_corr;
    }
  }
}
