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

#ifndef TOPPIC_MERGE_FEATURE_PRSM_READER_HPP_
#define TOPPIC_MERGE_FEATURE_PRSM_READER_HPP_

#include <fstream>

#include "merge/feature_prsm.hpp"

namespace toppic {

class FeaturePrsmReader {
 public:
  FeaturePrsmReader(const std::string &file_name);

  void close();

  FeaturePrsmPtr readOneFeature(); 

  FeaturePrsmPtrVec readAllFeatures();

 private:
  std::string file_name_;
  std::ifstream input_;
};

typedef std::shared_ptr<FeaturePrsmReader> FeaturePrsmReaderPtr;

}  // namespace toppic
#endif
