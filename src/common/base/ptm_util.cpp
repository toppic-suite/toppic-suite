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

#include <fstream>
#include <stdexcept>
#include <string>

#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "common/base/ptm_base.hpp"
#include "common/base/ptm_util.hpp"

namespace toppic {

namespace ptm_util {

struct DuplicatePtmError: public std::runtime_error {
  DuplicatePtmError(const std::string &msg) : std::runtime_error(msg) {}
};

PtmPtrVec readPtmTxt(const std::string &file_name) {
  std::ifstream infile(file_name.c_str());
  if (!infile.is_open()) {
    LOG_ERROR("Variable PTM file " << file_name <<  " cannot be opened!");
    throw std::runtime_error("Variable PTM file " + file_name + " cannot be opened!");
  }

  PtmPtrVec ptm_vec;

  std::string line;
  while (std::getline(infile, line)) {
    if (line.empty() || line[0] == '#') continue;
    line = str_util::rmComment(line);
    if (line == "") continue;
    try {
      std::vector<std::string> l = str_util::split(line, ",");
      if (l.size() != 5) throw std::runtime_error("Invalid format in line: " + line);
      if (l[2] == "*" && l[3] == "any") throw std::runtime_error("Invalid format in line: " + line);

      std::string name = l[0];
      double mass = std::stod(l[1]);
      int unimod_id = std::stoi(l[4]);
      PtmPtr ptm_ptr = PtmBase::getPtmPtr(
          std::make_shared<Ptm>(name, name, mass, unimod_id)); 

      //check if same name exists in the ptm_vec
      for (size_t i = 0; i < ptm_vec.size(); i++) {
        if (ptm_vec[i]->getName() == name) {
          throw DuplicatePtmError("Duplicate PTM name: " + name);
        }
      }
      ptm_vec.push_back(ptm_ptr);
    } catch (DuplicatePtmError& e) {
      LOG_ERROR("Errors in the Variable PTM file: " << file_name);
      LOG_ERROR("Errors in the line: " << line);
      throw std::runtime_error("Mod file cannot have two ptms with the same name: " + file_name);
    } catch (std::exception& e) {
      LOG_ERROR("Errors in the Variable PTM file: " << file_name);
      LOG_ERROR("Errors in the line: " << line);
      throw std::runtime_error("Errors in the Variable PTM file: " + file_name);
    }
  }
  infile.close();

  return ptm_vec;
}

}  // namespace ptm_util

}  // namespace toppic
