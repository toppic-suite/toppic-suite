//
// Created by abbash on 9/21/22.
//

#ifndef TOPPIC_GET_ENV_COLL_SCORE_HPP
#define TOPPIC_GET_ENV_COLL_SCORE_HPP

#include <cmath>
#include "env_util.hpp"
#include "fdeep/model.hpp"
#include "topfd/feature_detect/spectrum/peak_matrix.hpp"
#include "topfd/feature_detect/env_collection/env_collection.hpp"
#include "topfd/envcnn/env_cnn.hpp"
#include "topfd/feature_detect/test_output_functions/write_out_files.hpp"


namespace toppic {
namespace env_coll_score {
    double get_env_coll_score(fdeep::model& model, std::vector<double >& data);
};
}


#endif //TOPPIC_GET_ENV_COLL_SCORE_HPP
