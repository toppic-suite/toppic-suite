// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#include <string>

#include "base/string_util.hpp"
#include "base/ptm_base.hpp"
#include "base/ptm_util.hpp"

namespace prot {

namespace ptm_util {

PtmPtrVec readPtmTxt(const std::string &file_name) {
  std::ifstream infile(file_name.c_str());
  if (!infile.is_open()) {
    std::cerr << "Error: variable PTM file "
        << file_name <<  " can not be opened" << std::endl;
    exit(EXIT_FAILURE);
  }

  PtmPtrVec ptm_vec;

  std::string line;
  while (std::getline(infile, line)) {
    if (line[0] == '#') continue;
    line = string_util::rmComment(line);
    if (line == "") continue;
    try {
      std::vector<std::string> l = string_util::split(line, ',');
      if (l.size() != 5) throw line;

      if (l[2] == "*" && l[3] == "any") throw line;

      ptm_vec.push_back(PtmBase::getPtmPtr(std::make_shared<Ptm>(l[0], l[0], std::stod(l[1]), std::stoi(l[4]))));
    } catch (char const* e) {
      std::cerr << "Errors in the Variable PTM file: "
          << file_name << std::endl
          << "Please check the line" << std::endl
          << "\t" << e << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  infile.close();

  return ptm_vec;
}

}  // namespace ModUtil

}  // namespace prot
