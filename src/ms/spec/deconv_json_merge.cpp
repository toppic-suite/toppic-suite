//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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
#include "common/util/file_util.hpp"
#include "ms/spec/deconv_json_merge.hpp"

namespace toppic {

int DeconvJsonMerge::MAX_SPEC_NUM_PER_FILE = 1000000;

DeconvJsonMerge::DeconvJsonMerge(
    const std::vector<std::string> &spec_file_names,
    const std::string &output_file_name):
    spec_file_names_(spec_file_names),
    output_file_name_(output_file_name) {
    }

void DeconvJsonMerge::process() {
  std::vector<std::string> ms1_folder_names;
  std::vector<std::string> ms2_folder_names;
  for (size_t i = 0; i < spec_file_names_.size(); i++) { 
    std::string file_name = file_util::basename(spec_file_names_[i]);
    std::string base_name = file_name.substr(0, file_name.length() - 4);
    std::string ms1_folder = base_name + "_html" 
        + file_util::getFileSeparator()
        + "topfd" + file_util::getFileSeparator() + "ms1_json";
    ms1_folder_names.push_back(ms1_folder);
    std::string ms2_folder = base_name + "_html" 
        + file_util::getFileSeparator()
        + "topfd" + file_util::getFileSeparator() + "ms2_json";
    ms2_folder_names.push_back(ms2_folder);
  }
  
  std::string ms1_output_folder = output_file_name_ + "_html" 
      + file_util::getFileSeparator() 
      + "topfd" + file_util::getFileSeparator() + "ms1_json";
  std::string ms2_output_folder = output_file_name_ + "_html" 
      + file_util::getFileSeparator() 
      + "topfd" + file_util::getFileSeparator() + "ms2_json";

  mergeFiles(ms1_folder_names, ms1_output_folder, 
             MAX_SPEC_NUM_PER_FILE); 
  mergeFiles(ms2_folder_names, ms2_output_folder, 
             MAX_SPEC_NUM_PER_FILE); 
}

void DeconvJsonMerge::mergeFiles(const std::vector<std::string> &spec_folder_list,
                                 const std::string &output_folder, 
                                 int max_num_per_file) {
  if (file_util::exists(output_folder)) {
    file_util::delDir(output_folder);
  }
  file_util::createFolder(output_folder);
  for (size_t i = 0; i < spec_folder_list.size(); i++) {
    int id_base = max_num_per_file * i;
    file_util::copyJsonDir(spec_folder_list[i], output_folder, id_base);
  }
}

} /* namespace toppic */
