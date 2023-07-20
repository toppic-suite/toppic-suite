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

#include <cmath>

#include "ms/msmap/ms_map_peak.hpp"
#include "topfd/envcnn/onnx_env_cnn.hpp"
#include "topfd/ecscore/envelope/seed_envelope.hpp"
#include "topfd/ecscore/env_set/env_set.hpp"
#include "topfd/ecscore/env_set/env_set_util.hpp"
#include "topfd/ecscore/score/comp_env_cnn_score.hpp"

namespace toppic {

namespace comp_env_cnn_score {

int getIndex(double mz, double min_mz, double bin_size) {
  double mz_diff = mz - min_mz;
  int bin_idx = int(mz_diff / bin_size);
  return bin_idx;
}

MsMapPeakPtrVec  getIntvPeakList(MsMapPtr matrix_ptr, EnvSetPtr env_set_ptr,
                                 int spec_id) {
  EnvPeakPtrVec peak_list = env_set_ptr->getSeedPtr()->getPeakPtrList();
  double min_theo_peak = std::round(peak_list[0]->getPosition() * 1000.0) / 1000.0;
  double max_theo_peak = std::round(peak_list[peak_list.size() - 1]->getPosition() * 1000.0) / 1000.0;
  int start_idx = matrix_ptr->getColIndex(min_theo_peak - 0.1);
  int end_idx = matrix_ptr->getColIndex(max_theo_peak + 0.1);
  MsMapPeakPtrVec intv_peak_list;
  for (int peak_idx = start_idx; peak_idx <= end_idx; peak_idx++) {
    MsMapPeakPtrVec bin_peaks = matrix_ptr->getBinPeakList(spec_id, peak_idx);
    for (const auto &peak: bin_peaks)
      if (peak->getPosition() >= (min_theo_peak - 0.1) 
          && peak->getPosition() <= (max_theo_peak + 0.1))
        intv_peak_list.push_back(peak);
  }
  return intv_peak_list;
}

std::vector<std::vector<float>> getEnvcnnInputMatrix(MsMapPtr matrix_ptr,
                                                     EnvCollPtr coll_ptr) {
  std::vector<std::vector<float>> data_matrix = onnx_env_cnn::initInputMatrix(); 
  SeedEnvelopePtr seed_ptr = coll_ptr->getSeedPtr();
  EnvSetPtr env_set_ptr = coll_ptr->getSeedEnvSet();
  std::vector<double> theo_mz = env_set_ptr->getSeedMzList();
  std::vector<double> theo_inte = env_set_ptr->getSeedInteList();
  std::vector<double> exp_dist = env_set_util::getAggregateEnvelopeMz(env_set_ptr);
  std::vector<double> exp_dist_inte = env_set_util::getAggregateEnvelopeInte(env_set_ptr);
  double inte_ratio = env_set_util::calcInteRatio(theo_inte, exp_dist_inte);
  std::vector<double> scaled_theo_inte;
  EnvPeakPtrVec peak_list = env_set_ptr->getSeedPtr()->getPeakPtrList();
  for (const auto &peak: peak_list)
    scaled_theo_inte.push_back(inte_ratio * peak->getIntensity());
  double max_scalled_theo_inte = *std::max_element(scaled_theo_inte.begin(), scaled_theo_inte.end());
  int num_peaks = scaled_theo_inte.size();
  std::vector<double> normalize_theo_inte;
  std::vector<double> normalize_exp_inte;
  for (int i = 0; i < num_peaks; i++) {
    normalize_theo_inte.push_back(scaled_theo_inte[i] / max_scalled_theo_inte);
    normalize_exp_inte.push_back(exp_dist_inte[i] / max_scalled_theo_inte);
  }

  /// populate matrix
  double mass_tole = 0.02;
  double bin_size = matrix_ptr->getBinSize();
  num_peaks = theo_mz.size();
  for (int idx = 0; idx < num_peaks; idx++) {
    double exp_inte = normalize_exp_inte[idx];
    if (exp_inte == 0) continue;
    double mass_diff = std::abs(theo_mz[idx] - exp_dist[idx]);
    double md = mass_diff;
    if (abs(md) > mass_tole) {
      int r = std::rand() % 2;
      if (r == 0) {
        md = - mass_tole;
      }
      else {
        md = mass_tole;
      }
    }
    int p_idx = getIndex(std::round(theo_mz[idx] * 1000.0) / 1000.0, theo_mz[0], bin_size) + 10;
    if (p_idx >= 300) break;
    data_matrix[0][p_idx] = normalize_theo_inte[idx];
    data_matrix[1][p_idx] = exp_inte;
    data_matrix[2][p_idx] = md;
    data_matrix[3][p_idx] = normalize_theo_inte[idx] - exp_inte;
  }

  /// Add noise in EnvCNN matrix
  std::vector<float> noise_arr(300, 0.0);
  std::vector<std::vector<double>> noise_distribution_list;
  std::vector<std::vector<double>> noise_inte_distribution_list;
  for (int spec_id = coll_ptr->getStartSpecId(); spec_id <= coll_ptr->getEndSpecId(); spec_id++) {
    MsMapPeakPtrVec intv_peak_list = getIntvPeakList(matrix_ptr, env_set_ptr, spec_id);
    std::vector<double> t_noise_distribution_list;
    std::vector<double> t_noise_inte_distribution_list;
    for (const auto &elem: intv_peak_list) {
      t_noise_distribution_list.push_back(elem->getPosition());
      t_noise_inte_distribution_list.push_back(elem->getIntensity());
    }
    noise_distribution_list.push_back(t_noise_distribution_list);
    noise_inte_distribution_list.push_back(t_noise_inte_distribution_list);
  }

  int num_spec = noise_distribution_list.size();
  for (int idx = 0; idx < num_spec; idx++) {
    int num_peaks = noise_distribution_list[idx].size();
    for (int j_idx = 0; j_idx < num_peaks; j_idx++) {
      int p_idx =
        getIndex(std::round(noise_distribution_list[idx][j_idx] * 1000.0) / 1000.0, theo_mz[0], bin_size) + 10;
      if (p_idx >= 300) break;
      /// check if peak has been used as the data peak.
      bool matched_peak = false;
      // to accomodate +2 and -2 bins - reason tolerance of 0.02
      for (int shift = -2; shift <= 2; shift++) {
        int shifted_idx = p_idx + shift;
        if ((shifted_idx < 300) && (shifted_idx >= 0) && (data_matrix[0][shifted_idx] != 0)) {
          matched_peak = true;
        }
      }
      if (! matched_peak) {
        double max_scaled_theo_inte = *std::max_element(scaled_theo_inte.begin(), scaled_theo_inte.end());
        noise_arr[p_idx] = noise_arr[p_idx] + (noise_inte_distribution_list[idx][j_idx] / max_scaled_theo_inte);
      }
    }
  }
  /// add noise to the data
  num_peaks = noise_arr.size();
  for (int idx = 0; idx < num_peaks; idx++)
    data_matrix[1][idx] = data_matrix[1][idx] + noise_arr[idx];

  return data_matrix;
}


double compEnvcnnScore(MsMapPtr matrix_ptr, EnvCollPtr coll_ptr) {
  std::vector<std::vector<float>> matrix = getEnvcnnInputMatrix(matrix_ptr, coll_ptr); 
  std::vector<float> tensor;
  for (size_t i = 0; i < matrix.size(); i++) {
    tensor.insert(std::end(tensor), std::begin(matrix[i]), std::end(matrix[i]));
  }
  std::vector<double> results = onnx_env_cnn::predict(1, tensor);
  return results[0];
}

}
}
