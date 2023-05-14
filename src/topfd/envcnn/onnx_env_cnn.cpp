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

#include <assert.h>

#include <cstdlib>
#include <algorithm>
#include <iterator>
#include <limits>
#include <iostream>
#include <cmath>

#include "onnx/onnxruntime_cxx_api.h"

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "ms/spec/baseline_util.hpp"
#include "topfd/envcnn/onnx_env_cnn.hpp"

namespace toppic {

namespace onnx_env_cnn {

Ort::Env *env;
Ort::SessionOptions *session_options;
Ort::Session *session;
int batch_size = 1024;
std::vector<const char*> input_node_names = {"input"};
std::vector<const char*> output_node_names = {"output"};

void initModel(const std::string &dir_name, int thread_num) {
  std::string file_name = dir_name 
    + file_util::getFileSeparator() + "envcnn_models"
    + file_util::getFileSeparator() + "envcnn_two_block.onnx";

  env = new Ort::Env(ORT_LOGGING_LEVEL_WARNING, "env_cnn");
  session_options = new Ort::SessionOptions();
  session_options->SetIntraOpNumThreads(thread_num);

#ifdef _WIN32
  const wchar_t* model_path = L"envcnn_two_block.onnx";
#else
  const char* model_path = file_name.c_str();
#endif

  session = new Ort::Session(*env, model_path, *session_options);
}

MatchEnvPtr2D getBatchEnv(MatchEnvPtrVec &ori_envs) {
  MatchEnvPtr2D result;
  size_t start = 0;
  while (start < ori_envs.size()) {
    size_t end = start + batch_size;
    MatchEnvPtrVec one_batch;
    size_t idx = start;
    while (idx < end && idx < ori_envs.size()) {
      one_batch.push_back(ori_envs[idx]);
      idx++;
    }
    result.push_back(one_batch);
    start = start + batch_size;
  }
  return result;
}

void extractTheoPeakData(EnvelopePtr &theo_env, 
                         std::vector<double> &theo_mass,
                         std::vector<double> &theo_intes) {
  for (int i = 0; i < theo_env->getPeakNum(); i++) {
    theo_mass.push_back(theo_env->getMz(i));
    theo_intes.push_back(theo_env->getIntensity(i));
  }
}

void getTheoEnvData(MatchEnvPtr &ori_env, std::vector<double> &theo_mass, 
                    std::vector<double> &theo_intes,
                    double &max_inte, double &theo_mono_mz) {
  EnvelopePtr theo_env = ori_env->getTheoEnvPtr();
  theo_mono_mz= theo_env->getMonoMz();
  extractTheoPeakData(theo_env, theo_mass, theo_intes);
  max_inte= *std::max_element(std::begin(theo_intes), std::end(theo_intes));
  for (double & theo_inte : theo_intes) {
    theo_inte = theo_inte / max_inte;
  }
}

void getExpIntervaledPeakData(const PeakPtrVec &peak_list, double real_mono_mz,
                              std::vector<double> &theo_mass, 
                              std::vector<double> &peak_mass,
                              std::vector<double> &peak_intes) {
  double max_theo_mass = *std::max_element(std::begin(theo_mass), std::end(theo_mass));
  double min_theo_mass = *std::min_element(std::begin(theo_mass), std::end(theo_mass));
  for (const auto & k : peak_list) {
    double temp_mz = k->getPosition();
    if (temp_mz >= real_mono_mz - 5 && temp_mz <= real_mono_mz + 5) {
      if ((temp_mz >= min_theo_mass - 0.1) && (temp_mz <= max_theo_mass + 0.1)){
        peak_mass.push_back(temp_mz);
        peak_intes.push_back(k->getIntensity());
      }
    }
  }
}

void getExpEnvData(const PeakPtrVec &peak_list, MatchEnvPtr &ori_env, 
                   std::vector<double> &theo_mass,
                   std::vector<double> &peak_mass, 
                   std::vector<double> &peak_intes) {
  RealEnvPtr real_env = ori_env->getRealEnvPtr();
  double real_mono_mz = real_env->getMonoMz();
  getExpIntervaledPeakData(peak_list, real_mono_mz, theo_mass, peak_mass, peak_intes);
}

std::vector<std::vector<float>> initializeMatrix(int row, int col) {
  std::vector<std::vector<float>> matrix(row);
  for (int i = 0; i < row; i++) {
    matrix[i] = std::vector<float>(col);
    for (int j = 0; j < col; j++)
      matrix[i][j] = 0;
  }
  return matrix;
}

std::vector<std::vector<float>> initInputMatrix() {
  int row = 4;
  int col = 300;
  std::vector<std::vector<float>> matrix(row);
  for (int i = 0; i < row; i++) {
    matrix[i] = std::vector<float>(col);
    for (int j = 0; j < col; j++)
      matrix[i][j] = 0;
  }
  return matrix;
}


void extractFeature(const std::vector<double> &theo_masses, const std::vector<double> &theo_intes,
                    const std::vector<double> &peak_masses, const std::vector<double> &peak_intes,
                    double max_inte, double tolerance, size_t k, double &t_peak_mass,
                    double &t_peak_inte, double &exp_inte, double &md, double inte_diff) {
  t_peak_mass= theo_masses[k];
  t_peak_inte= theo_intes[k];
  exp_inte= 0;
  inte_diff= t_peak_inte;
  md= 0;
  int exp_peak_idx = -1;
  double mass_diff;
  double old_mass_diff = std::numeric_limits<double>::max();
  for (size_t exp_id = 0; exp_id < peak_masses.size(); exp_id++) {
    double exp_peak = peak_masses[exp_id];
    mass_diff = t_peak_mass - exp_peak;
    if (mass_diff <= tolerance && mass_diff >= -tolerance) {
      if (abs(mass_diff) < abs(old_mass_diff)){
        exp_peak_idx = exp_id;
        old_mass_diff = mass_diff;
      }
    }
  }
  mass_diff = -t_peak_mass;
  if (exp_peak_idx > -1) {
    double exp_peak = peak_masses[exp_peak_idx];
    exp_inte = peak_intes[exp_peak_idx] / max_inte;
    mass_diff = t_peak_mass - exp_peak;
    inte_diff = t_peak_inte - exp_inte;
  }
  // if mass difference is larger than 0.02, set it to
  // 0.02 or -0.02
  md = mass_diff;
  if (abs(md) > 0.02) {
    int r = std::rand() % 2;
    if (r == 0) {
      md = - 0.02;
    }
    else {
      md = 0.02;
    }
  }
}

void fillMatrix(std::vector<std::vector<float>> &matrix, double min_mz,
                double t_peak_mass, double t_peak_inte, double exp_inte,
                double mass_diff, double inte_diff) {
  int bin_index = int((t_peak_mass - min_mz) * 100);
  if (bin_index < 300 && bin_index >= 0 ){
    matrix[0][bin_index] = t_peak_inte;
    matrix[1][bin_index] = exp_inte;
    matrix[2][bin_index] = mass_diff;
    matrix[3][bin_index] = inte_diff;
  }
}

void addNoisePeaksToMatrix(std::vector<std::vector<float>> &matrix,
                           double min_mz, double max_inte,
                           const std::vector<double> &peak_masses,
                           const std::vector<double> &peak_intes) {
  for (size_t exp_id = 0; exp_id < peak_masses.size(); exp_id++) {
    double exp_peak = peak_masses[exp_id];
    int bin_index = int((exp_peak - min_mz) * 100);
    // Evaluate Peak Condition
    bool matched_peak = false;
    // to accomodate +2 and -2 bins - reason tolerance of 0.02
    for (int shift = -2; shift <= 2; shift++) {
      int shifted_idx = bin_index + shift;
      if ((shifted_idx < 300) && (shifted_idx >= 0) && (matrix[0][shifted_idx] != 0)) {
        matched_peak = true;
      }
    }
    if (!matched_peak && bin_index >= 0 && bin_index < 300) {
      matrix[1][bin_index] = (peak_intes[exp_id] / max_inte);
    }
  }
}

std::vector<float> getTensor(const PeakPtrVec &peak_list, MatchEnvPtr env) {
  std::vector<float> result;
  // extract theoretical mass and intensities into separate vectors
  std::vector<double> theo_masses;
  std::vector<double> theo_intes;
  double max_inte;
  double theo_mono_mz;
  getTheoEnvData(env, theo_masses, theo_intes, max_inte, theo_mono_mz);

  // extract intervaled experimental mass and intensities into separate vectors
  std::vector<double> peak_masses;
  std::vector<double> peak_intes;
  getExpEnvData(peak_list, env, theo_masses, peak_masses, peak_intes);

  // init matrix
  size_t row = 4; 
  size_t col = 300; 
  std::vector<std::vector<float>> matrix  = initializeMatrix(row, col);

  double tolerance = 0.02;

  double min_mz = theo_mono_mz - 0.1;
  for (size_t k = 0; k < theo_masses.size(); k++) {
    double t_peak_mass = 0;
    double t_peak_inte = 0;
    double exp_inte = 0;
    double inte_diff = 0;
    double mass_diff = 0;
    extractFeature(theo_masses, theo_intes, peak_masses, peak_intes, max_inte, tolerance, k, t_peak_mass,
                   t_peak_inte, exp_inte, mass_diff, inte_diff);
    fillMatrix(matrix, min_mz, t_peak_mass, t_peak_inte, exp_inte, mass_diff, inte_diff);
  }

  addNoisePeaksToMatrix(matrix, min_mz, max_inte, peak_masses, peak_intes);

  for (size_t i = 0; i < row; i++) {
    result.insert(std::end(result), std::begin(matrix[i]), std::end(matrix[i]));
  }
  return result;
}

std::vector<float> getBatchTensor(const PeakPtrVec &peak_list, MatchEnvPtrVec &envs){
  std::vector<float> result;
  for (size_t i = 0; i < envs.size(); i++) {
    std::vector<float> env_result = getTensor(peak_list, envs[i]); 
    result.insert(std::end(result), std::begin(env_result), std::end(env_result));
  }
  return result;
}

void predict(MatchEnvPtrVec &envs, std::vector<float> &input_tensor_values) {
  int env_num = envs.size();
  std::vector<double> pred_results = predict(env_num, input_tensor_values);

  for (int i = 0; i < env_num; i++) {
    envs[i]->setEnvcnnScore(pred_results[i]); 
    LOG_DEBUG(" msdeconv " << envs[i]->getMsdeconvScore() << " predict " << pred_results[i]);  
  }
}

std::vector<double> predict(int env_num, std::vector<float> &input_tensor_values) {
  std::vector<int64_t> input_node_dims {env_num, 4, 300};  
  size_t input_tensor_size; 
  input_tensor_size = env_num * 4 * 300;  

  // create input tensor object from data values
  auto memory_info = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);
  auto input_tensor = Ort::Value::CreateTensor<float>(memory_info,
                                                      input_tensor_values.data(),
                                                      input_tensor_size,
                                                      input_node_dims.data(), 3);
  assert(input_tensor.IsTensor());

  // score model & input tensor, get back output tensor
  auto output_tensors =
    session->Run(Ort::RunOptions{nullptr}, input_node_names.data(),
                 &input_tensor, 1, output_node_names.data(), 1);

  // Get pointer to output tensor float values
  float* float_arr = output_tensors.front().GetTensorMutableData<float>();

  std::vector<double> pred_results;
  for (int i = 0; i < env_num; i++) {
    double zero = float_arr[2*i];
    double one = float_arr[2*i+1];
    //softmax
    double score = std::exp(one)/ (std::exp(zero) + std::exp(one)); 
    pred_results.push_back(score);
  }
  return pred_results;
}

void computeEnvScores(PeakPtrVec &peak_list, MatchEnvPtrVec &ori_envs) {
  MatchEnvPtr2D batch_envs = getBatchEnv(ori_envs);

  for (size_t i = 0; i < batch_envs.size(); i++) {
    MatchEnvPtrVec envs = batch_envs[i];
    std::vector<float> tensor_data = getBatchTensor(peak_list, envs);
    predict(envs, tensor_data);
  }
}

}

}

