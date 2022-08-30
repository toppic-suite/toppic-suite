//
// Created by abbash on 8/26/22.
//

#ifndef TOPPIC_GET_ENV_CNN_SCORE_HPP
#define TOPPIC_GET_ENV_CNN_SCORE_HPP

#include "env_util.hpp"
#include "fdeep/model.hpp"
#include "topfd/feature_detect/spectrum/peak_matrix.hpp"
#include "topfd/feature_detect/env_collection/env_collection.hpp"
#include "topfd/envcnn/env_cnn.hpp"

namespace toppic {
namespace env_cnn_score {
    double get_envcnn_score(fdeep::model model, PeakMatrix peak_matrix, EnvCollection env_coll, double noiseIntensityLevel);
    std::vector<std::vector<double>> get_data_matrix_EnvCNN_aggregate_sum(PeakMatrix peak_matrix, EnvCollection env_coll, double noiseIntensityLevel, double bin_size);
}
}

#endif //TOPPIC_GET_ENV_CNN_SCORE_HPP
