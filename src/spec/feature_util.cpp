// Copyright (c) 2014 - 2018, The Trustees of Indiana University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>

#include "spec/feature_util.hpp"

namespace prot {

namespace feature_util {

void mergeFeatureFiles(const std::vector<std::string> & feature_file_lst,                                       
                       int N, const std::string & output_file) {
  std::ofstream outfile(output_file);
  for (size_t i = 0; i < feature_file_lst.size(); i++) {
    std::ifstream infile(feature_file_lst[i]);
    std::string line;
    while (std::getline(infile, line)) {
      if (line[0] == '#' || line == "" || line[0] == 'I') {
        outfile << line << std::endl;
        continue;
      }

      std::vector<std::string> strs;
      boost::split(strs, line, boost::is_any_of("\t "));
      outfile << N * i + std::stoi(strs[0]) << "\t";
      outfile << strs[1] << "\t";
      outfile << N * i + std::stoi(strs[2]) << "\t";
      outfile << strs[3] << "\t";
      for (size_t k = 4; k < strs.size(); k++) {
        outfile << strs[k] << "\t";
      }
      outfile << std::endl;
    }
  }
  outfile.close();
}

}  // namespace feature_util

}  // namespace prot
