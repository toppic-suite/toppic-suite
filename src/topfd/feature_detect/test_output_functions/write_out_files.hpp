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


#ifndef TOPPIC_WRITE_OUT_FILES_HPP
#define TOPPIC_WRITE_OUT_FILES_HPP

#include <fstream>
#include <string>
#include "topfd/feature_detect/spectrum/peak_matrix.hpp"
#include "topfd/feature_detect/envelope/seed_envelope.hpp"
#include "topfd/feature_detect/env_set/env_set.hpp"

namespace toppic {
  namespace write_out_files {
    void write_peak_matrix(PeakMatrix &peak_matrix, std::string file_name);

    void write_seed_envelopes(std::vector<SeedEnvelope> &seed_envs, std::string file_name);

    void write_noise_levels(PeakMatrix &peak_matrix, std::vector<double> &spec_noise_levels, std::string file_name);

    void write_env_set(PeakMatrix &peakMatrix, EnvSet &env_set, std::string file_name);

    void write_env_cnn_matrix(std::vector<std::vector<double>> &envcnn_data_matrix);

    void write_feature_map(std::vector<std::vector<double>> &matrix);

    void write_base_msz(std::vector<std::vector<double>> &matrix);
  }
}


#endif //TOPPIC_WRITE_OUT_FILES_HPP
