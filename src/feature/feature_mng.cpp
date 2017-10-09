// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <fstream>

#include <boost/algorithm/string.hpp>

#include "base/logger.hpp"
#include "feature/feature_mng.hpp" 

namespace prot {

EnvBasePtr FeatureMng::env_base_ptr_ = nullptr;

FeatureMng::FeatureMng(const std::string & exec_dir) {
  if (env_base_ptr_ == nullptr) {
    distr_file_name_ = exec_dir + distr_file_name_;
    LOG_DEBUG("distribution file name " << distr_file_name_);
    env_base_ptr_
        = std::make_shared<EnvBase>(distr_file_name_, distr_entry_num_, distr_mass_interval_);
    LOG_DEBUG("env base inited");
    env_rescore_para_file_name_ = exec_dir + env_rescore_para_file_name_;
    std::ifstream infile(env_rescore_para_file_name_);
    std::string line;
    while (std::getline(infile, line)) {
      std::vector<std::string> strs;
      boost::split(strs, line, boost::is_any_of("\t"));
      std::vector<double> scr;
      for (size_t i = 0; i < strs.size(); i++) {
        scr.push_back(std::stod(strs[i]));
      }
      env_rescore_para_.push_back(scr);
    }
    infile.close();
  }
}

void FeatureMng::setMinInte(double min_inte) {
  min_inte_ = min_inte;
  min_refer_inte_ = min_inte * sn_ratio_;
}

// get the mass group based on mass value 
int FeatureMng::getMassGroup(double base_mass) {
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
int FeatureMng::compMinConsPeakNum(int peak_num, int mass_group) {
  int min_cons_peak_num = peak_num + relative_consecutive_peak_num_;
  if (min_cons_peak_num < min_consecutive_peak_num_[mass_group]) {
    min_cons_peak_num = min_consecutive_peak_num_[mass_group];
  }
  return min_cons_peak_num;
}

void FeatureMng::setTolerance(double tolerance) {
  mz_tolerance_ = tolerance;
  score_error_tolerance_ = tolerance;
}

}
