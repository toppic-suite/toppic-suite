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

#include "common/util/str_util.hpp"
#include "quant/feature_prsm.hpp"

namespace toppic {

FeaturePrsm::FeaturePrsm(std::string line) {
  std::vector<std::string> strs;
  strs = str_util::split(line, "\t");
  id_ = std::stoi(strs[0]);
  mono_mass_ = std::stod(strs[1]);
  intensity_ = std::stod(strs[2]);
  retent_begin_ = std::stod(strs[3]);
  retent_end_ = std::stod(strs[4]);
  scan_begin_ = std::stoi(strs[5]);
  scan_end_ = std::stoi(strs[6]);
  min_charge_ = std::stoi(strs[7]);
  max_charge_ = std::stoi(strs[8]);
  if (strs.size() > 9) {
    prot_name_ = strs[9];
    prot_desc_ = strs[10];
    first_residue_ = std::stoi(strs[11]) - 1;
    last_residue_ = std::stoi(strs[12]) - 1;
    proteoform_ = strs[13];
    ms2_scan_ = std::stoi(strs[14]);
    prec_mass_ = std::stod(strs[15]);
  }
  else {
    ms2_scan_ = 0;
  }
}

}

