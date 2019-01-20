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

#ifndef TOPPIC_MERGE_FEATURE_SAMPLE_MERGE_HPP_
#define TOPPIC_MERGE_FEATURE_SAMPLE_MERGE_HPP_

#include "merge/feature_prsm.hpp"

namespace toppic {

class FeatureSampleMerge {
 public:
  FeatureSampleMerge(const std::vector<std::string> &input_file_names,
                     const std::string &output_file_name,
                     const std::string &db_file_name,
                     const std::string &fix_mod_str, 
                     double ppm);

  void process();

  void outputTable(FeaturePrsmPtrVec2D &table,
                   FeaturePrsmPtrVec &examples,
                   int sample_num);

 private:
  std::vector<std::string> input_file_names_;
  std::string output_file_name_;
  std::string db_file_name_;
  std::string fix_mod_str_;
  double ppm_;
};

typedef std::shared_ptr<FeatureSampleMerge> FeatureSampleMergePtr;

}  // namespace toppic

#endif 
