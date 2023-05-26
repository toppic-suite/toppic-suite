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
#include <cmath>

#include "common/util/logger.hpp"
#include "onnx/onnxruntime_cxx_api.h"

#include "common/util/file_util.hpp"
#include "topfd/ecscore/score/onnx_ecscore.hpp"


namespace toppic {

namespace onnx_ecscore {

Ort::Env *env;
Ort::SessionOptions *session_options;
Ort::Session *session;
int batch_size = 32;
std::vector<const char*> input_node_names = {"input"};
std::vector<const char*> output_node_names = {"output"};

void initModel(const std::string &dir_name, int thread_num) {
  std::string file_name = dir_name 
    + file_util::getFileSeparator() + "ecscore_models"
    + file_util::getFileSeparator() + "ecscore_seven_attr.onnx";

  env = new Ort::Env(ORT_LOGGING_LEVEL_WARNING, "ecsore");
  session_options = new Ort::SessionOptions();
  session_options->SetIntraOpNumThreads(thread_num);

#ifdef _WIN32
  const wchar_t* model_path = L"ecscore_seven_attr.onnx";
#else
  const char* model_path = file_name.c_str();
#endif

  session = new Ort::Session(*env, model_path, *session_options);
}


std::vector<double> predict(int feat_num, std::vector<float> &input_tensor_values) {
  std::vector<int64_t> input_node_dims {feat_num, 7};  
  size_t input_tensor_size; 
  input_tensor_size = feat_num * 7;  

  // create input tensor object from data values
  auto memory_info = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);
  auto input_tensor = Ort::Value::CreateTensor<float>(memory_info,
                                                      input_tensor_values.data(),
                                                      input_tensor_size,
                                                      input_node_dims.data(), 2);
  assert(input_tensor.IsTensor());

  // score model & input tensor, get back output tensor
  auto output_tensors =
    session->Run(Ort::RunOptions{nullptr}, input_node_names.data(),
                 &input_tensor, 1, output_node_names.data(), 1);

  // Get pointer to output tensor float values
  float* float_arr = output_tensors.front().GetTensorMutableData<float>();

  std::vector<double> pred_results;
  for (int i = 0; i < feat_num; i++) {
    double zero = float_arr[2*i];
    double one = float_arr[2*i+1];
    //softmax
    double score = std::exp(one)/ (std::exp(zero) + std::exp(one)); 
    pred_results.push_back(score);
  }
  return pred_results;
}

double predict(std::vector<float> &input_tensor_values) {
  std::vector<double> results = predict(1, input_tensor_values);
  return results[0];
}

}
}
