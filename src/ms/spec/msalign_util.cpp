//Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane University.
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

#include "ms/spec/msalign_util.hpp"

#include <fstream>
#include <string>

#include "common/util/logger.hpp"
#include "ms/spec/msalign_reader.hpp"

namespace toppic {

namespace msalign_util {

int countSpNum(const std::string &spectrum_file_name) {
  MsAlignReader reader(spectrum_file_name);
  int cnt = 0;
  DeconvMsPtr deconv_ms_ptr;
  while ((deconv_ms_ptr = reader.getNextMsPtr()) != nullptr) {
    cnt++;
  }
  return cnt;
}

void geneSpIndex(const std::string &spectrum_file_name) {
  int sp_num = countSpNum(spectrum_file_name);
  std::string index_file_name = spectrum_file_name + "_index";
  std::ofstream index_output(index_file_name);
  index_output << sp_num << std::endl;
}

int getSpNum(const std::string &spectrum_file_name) {
  std::string index_file_name = spectrum_file_name + "_index";
  std::ifstream index_input(index_file_name);
  std::string line;
  std::getline(index_input, line);
  int sp_num = std::stoi(line);
  LOG_DEBUG("Get sp number " << sp_num);
  return sp_num;
}

}  // namespace msalign_util

}  // namespace toppic
