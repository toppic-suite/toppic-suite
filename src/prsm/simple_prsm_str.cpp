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

#include <string>
#include <vector>

#include "common/util/logger.hpp"
#include "prsm/prsm_util.hpp"
#include "prsm/simple_prsm_str.hpp"

namespace toppic {

SimplePrsmStr::SimplePrsmStr(const std::vector<std::string> &str_vec) {
  str_vec_ = str_vec;
  std::string line = prsm_util::getXmlLine(str_vec_, "<file_name>");
  file_name_ = prsm_util::getValueStr(line);
  line = prsm_util::getXmlLine(str_vec_, "<spectrum_id>");
  spectrum_id_ = std::stoi(prsm_util::getValueStr(line));
  line = prsm_util::getXmlLine(str_vec_, "<score>");
  score_ = std::stod(prsm_util::getValueStr(line));
  line = prsm_util::getXmlLine(str_vec_, "<sequence_name>");
  seq_name_ = prsm_util::getValueStr(line);
  line = prsm_util::getXmlLine(str_vec_, "<sequence_desc>");
  seq_desc_ = prsm_util::getValueStr(line);
}

}  // namespace toppic
