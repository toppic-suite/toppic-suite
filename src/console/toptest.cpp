#include <iostream>

#include <assert.h>
#include "onnx/onnxruntime_cxx_api.h"

int main(int argc, char* argv[]) {

  Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "env_cnn");

  Ort::SessionOptions session_options;
  session_options.SetIntraOpNumThreads(1);

#ifdef _WIN32
  const wchar_t* model_path = L"envcnn_two_block.onnx";
#else
  const char* model_path = "envcnn_two_block.onnx";
#endif

  Ort::Session session(env, model_path, session_options);

  /**
  // print model input layer (node names, types, shape etc.)
  Ort::AllocatorWithDefaultOptions allocator;

  // print number of model input nodes
  const size_t num_input_nodes = session.GetInputCount();
  std::vector<Ort::AllocatedStringPtr> input_names_ptr;
  std::vector<const char*> input_node_names;
  input_names_ptr.reserve(num_input_nodes);
  input_node_names.reserve(num_input_nodes);
  std::vector<int64_t> input_node_dims;  // simplify... this model has only 1
                                         // input node {1, 4, 300}.
                                         // Otherwise need vector<vector<>>

  std::cout << "Number of inputs = " << num_input_nodes << std::endl;

  // iterate over all input nodes
  for (size_t i = 0; i < num_input_nodes; i++) {
    // print input node names
    auto input_name = session.GetInputNameAllocated(i, allocator);
    std::cout << "Input " << i << " : name =" << input_name.get() << std::endl;
    input_node_names.push_back(input_name.get());
    input_names_ptr.push_back(std::move(input_name));

    // print input node types
    auto type_info = session.GetInputTypeInfo(i);
    auto tensor_info = type_info.GetTensorTypeAndShapeInfo();

    ONNXTensorElementDataType type = tensor_info.GetElementType();
    std::cout << "Input " << i << " : type = " << type << std::endl;

    // print input shapes/dims
    input_node_dims = tensor_info.GetShape();
    std::cout << "Input " << i << " : num_dims = " << input_node_dims.size() << '\n';
    for (size_t j = 0; j < input_node_dims.size(); j++) {
      std::cout << "Input " << i << " : dim[" << j << "] =" << input_node_dims[j] << '\n';
    }
    std::cout << std::flush;
  }
  */

  std::vector<const char*> input_node_names = {"input"};
  std::vector<int64_t> input_node_dims {2, 4, 300};  
  constexpr size_t input_tensor_size = 2 * 4 * 300;  

  std::vector<float> input_tensor_values(input_tensor_size);
  std::vector<const char*> output_node_names = {"output"};

  // initialize input data with values in [0.0, 1.0]
  for (unsigned int i = 0; i < input_tensor_size/2; i++) {
    input_tensor_values[i] = (float)i / (input_tensor_size + 1);
    input_tensor_values[input_tensor_size/2 + i] = (float)i / (input_tensor_size + 1);
  }

  // create input tensor object from data values
  auto memory_info = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);
  auto input_tensor = Ort::Value::CreateTensor<float>(memory_info, input_tensor_values.data(), input_tensor_size,
                                                            input_node_dims.data(), 3);
  assert(input_tensor.IsTensor());

  // score model & input tensor, get back output tensor
  auto output_tensors =
      session.Run(Ort::RunOptions{nullptr}, input_node_names.data(),
                  &input_tensor, 1, output_node_names.data(), 1);

  // Get pointer to output tensor float values
  float* floatarr = output_tensors.front().GetTensorMutableData<float>();

  // score the model, and print scores for first 5 classes
  for (int i = 0; i < 4; i++) {
    std::cout << "Score for class [" << i << "] =  " << floatarr[i] << '\n';
  }
  std::cout << std::flush;

  std::cout << "Done!" << std::endl;

  return 0;
}

