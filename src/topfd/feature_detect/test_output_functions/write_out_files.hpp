//
// Created by abbash on 9/7/22.
//

#ifndef TOPPIC_WRITE_OUT_FILES_HPP
#define TOPPIC_WRITE_OUT_FILES_HPP

#include <fstream>
#include <string>
#include "topfd/feature_detect/spectrum/peak_matrix.hpp"
#include "topfd/feature_detect/envelope/seed_envelope.hpp"
#include "topfd/feature_detect/env_set/env_set.hpp"

namespace toppic {
namespace write_out_files {
    void write_peak_matrix(PeakMatrix& peak_matrix, std::string file_name);
    void write_seed_envelopes(std::vector<SeedEnvelope>& seed_envs, std::string file_name);
    void write_noise_levels(PeakMatrix& peak_matrix, std::vector<double>& spec_noise_levels, std::string file_name);
    void write_env_set(PeakMatrix& peakMatrix, EnvSet& env_set, std::string file_name);
    void write_env_cnn_matrix(std::vector<std::vector<double>>& envcnn_data_matrix);
}
}



#endif //TOPPIC_WRITE_OUT_FILES_HPP
