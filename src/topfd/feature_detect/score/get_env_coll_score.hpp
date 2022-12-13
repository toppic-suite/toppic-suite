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
    double get_env_coll_score(fdeep::model &model, std::vector<double> &data);
  };
}

#endif //TOPPIC_GET_ENV_COLL_SCORE_HPP
