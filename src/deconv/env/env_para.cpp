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

#include "util/logger.hpp"
#include "deconv/env/env_para.hpp" 

namespace toppic {

void EnvPara::setMinInte(double min_inte, int ms_level) {
  min_inte_ = min_inte;
  if (ms_level == 1) {
    min_refer_inte_ = min_inte * ms_one_sn_ratio_;
  }
  else {
    min_refer_inte_ = min_inte * ms_two_sn_ratio_;
  }
}

// get the mass group based on mass value 
int EnvPara::getMassGroup(double base_mass) {
  int group = -1;
  for (int i = 0; i < (int)mass_group_boundary_.size() - 1; i++) {
    if (base_mass >= mass_group_boundary_[i] && base_mass < mass_group_boundary_[i + 1]) {
      group = i;
      break;
    }
  }
  return group;
}

// compute minimum consecutive peak num: max{peak_num -3 , 3} 
int EnvPara::compMinConsPeakNum(int peak_num, int mass_group) {
  int min_cons_peak_num = peak_num + relative_consecutive_peak_num_;
  if (min_cons_peak_num < min_consecutive_peak_num_[mass_group]) {
    min_cons_peak_num = min_consecutive_peak_num_[mass_group];
  }
  return min_cons_peak_num;
}

void EnvPara::setTolerance(double tolerance) {
  mz_tolerance_ = tolerance;
  score_error_tolerance_ = tolerance; 
}

  /*  
  env_rescore_para_file_name_ = resource_dir + env_rescore_para_file_name_;
  std::ifstream infile(env_rescore_para_file_name_);
  std::string line;
  while (std::getline(infile, line)) {
    //boost::split(strs, line, boost::is_any_of("\t"));
    std::vector<std::string> strs = str_util::split(line, "\t");
    std::vector<double> scr;
    for (size_t i = 0; i < strs.size(); i++) {
      scr.push_back(std::stod(strs[i]));
    }
    env_rescore_para_.push_back(scr);
  }
  infile.close();
  */

}
