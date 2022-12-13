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

#include "get_env_cnn_score.hpp"

namespace toppic {
  namespace env_cnn_score {
    int get_index(double mz, double min_mz, double bin_size) {
      double mz_diff = mz - min_mz;
      int bin_idx = int(mz_diff / bin_size);
      return bin_idx;
    }

    double round(double val) {
      if (val < 0)
        return ceil(val - 0.5);
      return floor(val + 0.5);
    }

    std::vector<ExpPeak> get_intv_peak_list(PeakMatrix &peak_matrix, EnvSet &env_set, int spec_id) {
      std::vector<SimplePeak> peak_list = env_set.get_peak_list();
      double min_theo_peak = round(peak_list[0].getPos() * 1000.0) / 1000.0;
      double max_theo_peak = round(peak_list[peak_list.size() - 1].getPos() * 1000.0) / 1000.0;
      int start_idx = peak_matrix.get_index(min_theo_peak - 0.1);
      int end_idx = peak_matrix.get_index(max_theo_peak + 0.1);
      std::vector<ExpPeak> intv_peak_list;
      for (int peak_idx = start_idx; peak_idx <= end_idx; peak_idx++) {
        std::vector<ExpPeak> bin_peaks = peak_matrix.get_bin_peak(spec_id, peak_idx);
        for (const auto &peak: bin_peaks)
          if (peak.getPos() >= (min_theo_peak - 0.1) && peak.getPos() <= (max_theo_peak + 0.1))
            intv_peak_list.push_back(peak);
      }
      return intv_peak_list;
    }

    double get_envcnn_score(fdeep::model &model, PeakMatrix &peak_matrix, EnvCollection &env_coll,
                            double noiseIntensityLevel) {
      std::vector<std::vector<double>> envcnn_data_matrix = get_data_matrix_EnvCNN_aggregate_sum(peak_matrix, env_coll,
                                                                                                 noiseIntensityLevel,
                                                                                                 0.01);
      fdeep::tensor_shape tensor_shape(300, 5);
      fdeep::tensor t(tensor_shape, 0.0f);
      for (int y = 0; y < 300; ++y)
        for (int x = 0; x < 5; ++x)
          t.set(0, 0, 0, y, x, envcnn_data_matrix[y][x]);
      const auto result = model.predict({t});
      return result[0].get(0, 0, 0, 0, 0);
    }

    std::vector<std::vector<double>>
    get_data_matrix_EnvCNN_aggregate_sum(PeakMatrix &peak_matrix, EnvCollection &env_coll,
                                         double noiseIntensityLevel, double bin_size) {
      SeedEnvelope seed_env = env_coll.getSeedEnv();
      EnvSet env_set = env_coll.get_seed_env_set();
      std::vector<double> theo_mz = env_set.get_theo_distribution_mz();
      std::vector<double> theo_inte = env_set.get_theo_distribution_inte();
      double mass_tole;
      std::vector<std::vector<double>> envcnn_data_matrix = env_cnn::initializeMatrix(mass_tole);
      std::vector<double> exp_dist = env_utils::get_aggregate_envelopes_mz(env_set);
      std::vector<double> exp_dist_inte = env_utils::get_aggregate_envelopes_inte(env_set);
      double inte_ratio = env_utils::calcInteRatio_scan(theo_inte, exp_dist_inte);
      std::vector<double> scaled_theo_inte;
      std::vector<SimplePeak> peak_list = env_set.get_peak_list();
      for (const auto &peak: peak_list)
        scaled_theo_inte.push_back(inte_ratio * peak.getInte());
      double max_scalled_theo_inte = *std::max_element(scaled_theo_inte.begin(), scaled_theo_inte.end());

      int num_peaks = scaled_theo_inte.size();
      std::vector<double> normalize_theo_inte;
      std::vector<double> normalize_exp_inte;
      for (int i = 0; i < num_peaks; i++) {
        normalize_theo_inte.push_back(scaled_theo_inte[i] / max_scalled_theo_inte);
        normalize_exp_inte.push_back(exp_dist_inte[i] / max_scalled_theo_inte);
      }

      /// populate matrix
      num_peaks = theo_mz.size();
      for (int idx = 0; idx < num_peaks; idx++) {
        double exp_inte = normalize_exp_inte[idx];
        if (exp_inte == 0) continue;
        double mass_diff = std::abs(theo_mz[idx] - exp_dist[idx]);
        double mass_diff_score = 0;
        if (mass_diff < mass_tole)
          mass_diff_score = (mass_tole - mass_diff) / mass_tole;
        int p_idx = get_index(round(theo_mz[idx] * 1000.0) / 1000.0, theo_mz[0], bin_size) + 10;
        if (p_idx >= 300) break;
        envcnn_data_matrix[p_idx][0] = normalize_theo_inte[idx];
        envcnn_data_matrix[p_idx][1] = exp_inte;
        envcnn_data_matrix[p_idx][2] = mass_diff_score;
        envcnn_data_matrix[p_idx][3] = normalize_theo_inte[idx] - exp_inte;
        envcnn_data_matrix[p_idx][4] = log10(
            *std::max_element(scaled_theo_inte.begin(), scaled_theo_inte.end()) / noiseIntensityLevel);
      }

      /// Add noise in EnvCNN matrix
      std::vector<double> noise_arr(300, 0.0);
      std::vector<std::vector<double>> noise_distribution_list;
      std::vector<std::vector<double>> noise_inte_distribution_list;
      for (int spec_id = env_coll.getStartSpecId(); spec_id <= env_coll.getEndSpecId(); spec_id++) {
        std::vector<ExpPeak> intv_peak_list = get_intv_peak_list(peak_matrix, env_set, spec_id);
        std::vector<double> t_noise_distribution_list;
        std::vector<double> t_noise_inte_distribution_list;
        for (const auto &elem: intv_peak_list) {
          t_noise_distribution_list.push_back(elem.getPos());
          t_noise_inte_distribution_list.push_back(elem.getInte());
        }
        noise_distribution_list.push_back(t_noise_distribution_list);
        noise_inte_distribution_list.push_back(t_noise_inte_distribution_list);
      }

      int num_spec = noise_distribution_list.size();
      for (int idx = 0; idx < num_spec; idx++) {
        int num_peaks = noise_distribution_list[idx].size();
        for (int j_idx = 0; j_idx < num_peaks; j_idx++) {
          int p_idx =
              get_index(round(noise_distribution_list[idx][j_idx] * 1000.0) / 1000.0, theo_mz[0], bin_size) + 10;
          if (p_idx >= 300) break;
          /// check if peak has been used as the data peak.
          bool peak_condition = false;
          for (int i = 0; i < 3; i++) {
            if ((p_idx + i < 300) && (p_idx - i >= 0) && (envcnn_data_matrix[p_idx][0] == 0)) {
              if (envcnn_data_matrix[p_idx - i][0] == 0 and envcnn_data_matrix[p_idx + i][0] == 0)
                peak_condition = true;
              else {
                peak_condition = false;
                break;
              }
            }
          }
          if (peak_condition) {
            double max_scaled_theo_inte = *std::max_element(scaled_theo_inte.begin(), scaled_theo_inte.end());
            noise_arr[p_idx] = noise_arr[p_idx] + (noise_inte_distribution_list[idx][j_idx] / max_scaled_theo_inte);
          }
        }
      }
      /// add noise to the data
      num_peaks = noise_arr.size();
      for (int idx = 0; idx < num_peaks; idx++)
        envcnn_data_matrix[idx][1] = envcnn_data_matrix[idx][1] + noise_arr[idx];
      return envcnn_data_matrix;
    }
  }
}
