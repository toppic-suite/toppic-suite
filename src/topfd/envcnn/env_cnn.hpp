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

#ifndef TOPPIC_TOPFD_ENVCNN_ENV_CNN_HPP_
#define TOPPIC_TOPFD_ENVCNN_ENV_CNN_HPP_

#include <string>
#include "fdeep/fdeep.hpp"
#include "ms/env/match_env.hpp"

namespace toppic {

namespace env_cnn {

  void initModel(const std::string &dir_name, int thread_num);
  std::vector<std::vector<double>> initializeMatrix(double &tolerance);
  void generateTensors(std::vector<fdeep::tensors> &tensorsL, const std::vector<std::vector<double>> &matrix);
  void computeEnvScores(MatchEnvPtrVec &ori_env, PeakPtrVec &peak_list); 

}

}

#endif
