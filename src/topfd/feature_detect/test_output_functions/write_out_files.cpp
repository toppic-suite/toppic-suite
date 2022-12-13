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

#include "write_out_files.hpp"

namespace toppic {
  void write_out_files::write_peak_matrix(PeakMatrix &peak_matrix, std::string file_name) {
    std::ofstream out_file;
    int spec_num = peak_matrix.get_spec_num();
    for (int spec_id = 0; spec_id < spec_num; spec_id++) {
      out_file.open(file_name + "_" + std::to_string(spec_id));
      std::vector<std::vector<ExpPeak>> peak_row = peak_matrix.getRow(spec_id);
      int bin_arr_num = peak_row.size();
      for (int bin_arr_id = 0; bin_arr_id < bin_arr_num; bin_arr_id++) {
        std::vector<ExpPeak> bin_arr = peak_row[bin_arr_id];
        int num_bins = bin_arr.size();
        for (int bin_id = 0; bin_id < num_bins; bin_id++)
          if (!bin_arr[bin_id].isEmpty())
            out_file << "ID: (" << spec_id << ", " << bin_arr_id << ", " << bin_id << ") "
                     << bin_arr[bin_id].getString();
      }
      out_file.close();
    }
  }

  void write_out_files::write_seed_envelopes(std::vector<SeedEnvelope> &seed_envs, std::string file_name) {
    int num_seed_envs = seed_envs.size();
    std::ofstream out_file;
    out_file.open(file_name);
    for (int seed_env_id = 0; seed_env_id < num_seed_envs; seed_env_id++)
      out_file << "ID: " << seed_env_id << ", " << seed_envs[seed_env_id].getString();
    out_file.close();
  }

  void write_out_files::write_noise_levels(PeakMatrix &peak_matrix, std::vector<double> &spec_noise_levels,
                                           std::string file_name) {
    int num_spec = spec_noise_levels.size();
    std::ofstream out_file;
    out_file.open(file_name);
    out_file << "Data Level Noise intensity: " << peak_matrix.get_min_inte() << "\n";
    for (int id = 0; id < num_spec; id++)
      out_file << "Spectrum " << id << ": " << spec_noise_levels[id] << "\n";
    out_file.close();
  }

  void write_out_files::write_env_set(PeakMatrix &peakMatrix, EnvSet &env_set, std::string file_name) {
    double snr = 3.0;
    double noise_inte = peakMatrix.get_min_inte();

    std::ofstream out_file;
    out_file.open(file_name, std::ios_base::app);

    SeedEnvelope seed_env = env_set.getSeedEnv();
    std::vector<std::vector<double>> map = env_set.get_map(3.0, peakMatrix.get_min_inte());
    Xic xic = env_set.getXic();
    std::vector<ExpEnvelope> env_list = env_set.getExpEnvList();

    out_file << "SEED ENVELOPE \n";
    out_file << "SPEC ID: " << seed_env.getSpecId() << ", ";
    out_file << "Env ID: " << seed_env.getEnvId() << ", ";
    out_file << "Charge: " << seed_env.getCharge() << ", ";
    out_file << "Mass: " << seed_env.getMass() << ", ";
    out_file << "Pos: " << seed_env.getPos() << ", ";
    out_file << "Inte: " << env_set.comp_intensity(snr, noise_inte) << "\n";

    out_file << "Distribution: ";
    for (auto p: env_set.get_theo_distribution_mz())
      out_file << p << " ";
    out_file << "\n";
    out_file << "Distribution Inte: ";
    for (auto p: env_set.get_theo_distribution_inte())
      out_file << p << " ";
    out_file << "\n";
    out_file << "Feature Boundary: " << env_set.getStartSpecId() << ", " << env_set.getEndSpecId() << "\n";

    out_file << "Envelope XIC: ";
    for (auto s: xic.getInteList())
      out_file << s << " ";
    out_file << "\n";

    out_file << "Smoothed Envelope XIC: ";
    for (auto s: xic.getSmoothedInteList())
      out_file << s << " ";
    out_file << "\n";

    out_file << "Theo Map\n";
    for (auto s: map) {
      for (auto p: s)
        out_file << p << " ";
      out_file << "\n";
    }
    out_file << "\n\n";

    out_file << "Exp Map\n";
    for (auto env: env_list) {
      for (auto p: env.getExpEnvList())
        out_file << p.getInte() << " ";
      out_file << "\n";
    }
    out_file << "\n\n";

    out_file << "Exp Map - mz \n";
    for (auto env: env_list) {
      for (auto p: env.getExpEnvList())
        out_file << p.getPos() << " ";
      out_file << "\n";
    }
    out_file << "\n\n";

    out_file.close();
  }

  void write_out_files::write_env_cnn_matrix(std::vector<std::vector<double>> &envcnn_data_matrix) {
    std::ofstream out_file;
    out_file.open("envcnn.txt", std::ios_base::app);
    for (auto &row: envcnn_data_matrix) {
      for (auto &col: row)
        out_file << col << ", ";
      out_file << "\n";
    }
    out_file << "\n\n";
    out_file.close();
  }

  void write_out_files::write_feature_map(std::vector<std::vector<double>> &matrix) {
    std::ofstream out_file;
    out_file.open("matrix.csv", std::ios_base::app);
    for (auto &row: matrix) {
      for (auto &col: row)
        out_file << col << ", ";
      out_file << "\n";
    }
    out_file << "\n\n";
    out_file.close();
  }

  void write_out_files::write_base_msz(std::vector<std::vector<double>> &matrix) {
    std::ofstream out_file;
    out_file.open("base.csv", std::ios_base::app);
    for (auto &row: matrix) {
      for (auto &col: row)
        out_file << col << ", ";
      out_file << "\n";
    }
    out_file << "\n\n";
    out_file.close();
  }
}