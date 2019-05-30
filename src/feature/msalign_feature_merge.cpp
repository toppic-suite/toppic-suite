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

#include <set>
#include <algorithm>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "common/util/str_util.hpp"
#include "spec/msalign_frac_merge.hpp"
#include "feature/frac_feature_merge.hpp"
#include "feature/msalign_feature_merge.hpp"

namespace toppic {

int MsAlignFeatureMerge::MAX_NUM_PER_FILE = 1000000;

MsAlignFeatureMerge::MsAlignFeatureMerge(
    const std::vector<std::string> &spec_file_names,
    const std::string &output_file_name):
    spec_file_names_(spec_file_names),
    output_file_name_(output_file_name) {
    }

void MsAlignFeatureMerge::process(std::string &para_str) {
  std::vector<std::string> ms1_file_names;
  std::vector<std::string> ms2_file_names;
  std::vector<std::string> ms1_feature_names;
  std::vector<std::string> ms2_feature_names;
  for (size_t i = 0; i < spec_file_names_.size(); i++) { 
    std::string base_name = file_util::basename(spec_file_names_[i]);
    std::string ms1_name = base_name + "_ms1.msalign";
    ms1_file_names.push_back(ms1_name);
    std::string ms2_name = base_name + "_ms2.msalign";
    ms2_file_names.push_back(ms2_name);
    std::string ms1_feature = base_name + "_ms1.feature";
    ms1_feature_names.push_back(ms1_feature);
    std::string ms2_feature = base_name + "_ms2.feature";
    ms2_feature_names.push_back(ms2_feature);
  }
  
  std::string ms1_spec_output_name = output_file_name_ + "_ms1.msalign";
  std::string ms2_spec_output_name = output_file_name_ + "_ms2.msalign";
  std::string ms1_feature_output_name = output_file_name_ + "_ms1.feature";
  std::string ms2_feature_output_name = output_file_name_ + "_ms2.feature";

  msalign_frac_merge::mergeFiles(ms1_file_names, ms1_spec_output_name, 
                                 MAX_NUM_PER_FILE, para_str); 
  msalign_frac_merge::mergeFiles(ms2_file_names, ms2_spec_output_name, 
                                 MAX_NUM_PER_FILE, para_str); 
  frac_feature_merge::mergeFiles(ms1_feature_names, ms1_feature_output_name, 
                                 ms2_feature_names, ms2_feature_output_name,
                                 MAX_NUM_PER_FILE, para_str); 
}

} /* namespace toppic */
