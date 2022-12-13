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

#include "get_env_coll_score.hpp"

namespace toppic {
  void generateTensors(std::vector<fdeep::tensors> &tensorsL, const std::vector<double> &matrix) {
    fdeep::tensor_shape tensor_shape(8);
    fdeep::tensor t(tensor_shape, 0.0f);
    for (int y = 0; y < 8; ++y)
      t.set(0, 0, 0, 0, y, matrix[y]);
    std::vector<fdeep::tensor> tensors;
    tensors.push_back(t);
    tensorsL.push_back(tensors);
  }

  double env_coll_score::get_env_coll_score(fdeep::model &model, std::vector<double> &data) {
    std::vector<fdeep::tensors> tensorsL;
    generateTensors(tensorsL, data);
    if (!tensorsL.empty()) {
      std::vector<fdeep::tensors> pred_scores = model.predict_multi(tensorsL, false);
      return pred_scores[0][0].get(0, 0, 0, 0, 0);
    }
    return 0.0;
  }
}