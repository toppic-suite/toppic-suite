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

#ifndef TOPPIC_FEATURE_SPEC_FEATURE_WRITER_HPP_
#define TOPPIC_FEATURE_SPEC_FEATURE_WRITER_HPP_

#include <memory>
#include <vector>
#include <string>
#include <fstream>

#include "feature/spec_feature.hpp"

namespace toppic {

namespace spec_feature_writer {

void writeFeatures(const std::string &output_file_name,
                   const SpecFeaturePtrVec &features);

void writeHeader(std::ofstream &of); 

void writeOneFeature(std::ofstream &of, SpecFeaturePtr feature);

}

} /* namespace toppic */

#endif
