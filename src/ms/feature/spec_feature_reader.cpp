//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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

#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "ms/feature/spec_feature_reader.hpp"

namespace toppic {

SpecFeatureReader::SpecFeatureReader(const std::string &file_name):
    file_name_(file_name) {
      input_.open(file_name.c_str(), std::ios::in);
      if (!input_.is_open()) {
        LOG_ERROR("Feature file  " << file_name << " does not exist.");
        exit(EXIT_FAILURE); 
      }
      // read header line
      std::string line;
      std::getline(input_, line);
    }

SpecFeatureReader::~SpecFeatureReader() {
  if (input_.is_open()) {
    input_.close();
  }
}

void SpecFeatureReader::close() {
  input_.close();
}

SpecFeaturePtr SpecFeatureReader::readOneFeature() {
  std::string line; 
  if (std::getline(input_, line)) {
    str_util::trim(line);
    //std::cout << "line " << line << std::endl;
    SpecFeaturePtr feature = std::make_shared<SpecFeature>(line);
    //std::cout << "ms/feature created " << std::endl;
    return feature;
  }
  else {
    return nullptr;
  }
}


SpecFeaturePtrVec SpecFeatureReader::readAllFeatures() {
  SpecFeaturePtrVec all_features;
  SpecFeaturePtr feature;
  while ((feature = readOneFeature()) != nullptr) {
    all_features.push_back(feature);
  }
  return all_features;
}


}  // namespace prot
