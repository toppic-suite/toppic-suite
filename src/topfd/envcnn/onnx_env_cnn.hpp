//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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

#ifndef TOPPIC_TOPFD_ONNX_ENVCNN_ENV_CNN_HPP_
#define TOPPIC_TOPFD_ONNX_ENVCNN_ENV_CNN_HPP_

#include <string>
#include <vector>

#include "ms/env/match_env.hpp"

namespace toppic {

namespace onnx_env_cnn {

  void initModel(const std::string &dir_name, int thread_num);

  void computeEnvScores(PeakPtrVec &peak_list, MatchEnvPtrVec &ori_env); 

  std::vector<std::vector<float>> initInputMatrix();

  std::vector<double> predict(int env_num, std::vector<float> &input_tensor_values); 
}

}

#endif