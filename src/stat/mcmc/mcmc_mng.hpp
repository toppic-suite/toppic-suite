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


#ifndef TOPPIC_STAT_MCMC_MNG_HPP_
#define TOPPIC_STAT_MCMC_MNG_HPP_

#include <cmath>

#include "para/prsm_para.hpp"

namespace toppic {

class MCMCMng {
 public:
  MCMCMng(PrsmParaPtr prsm_para_ptr, 
          const std::string & input_file_ext, 
          const std::string & output_file_ext,
          const std::string & residue_mod_file,
          int max_known_mods,
          int thread_num):
      prsm_para_ptr_(prsm_para_ptr),
      input_file_ext_(input_file_ext),
      output_file_ext_(output_file_ext),
      residue_mod_file_(residue_mod_file),
      max_known_mods_(max_known_mods),
      thread_num_(thread_num) {
        convert_ratio_ = 274.335215;
        error_tolerance_ = 0.1;
        int_tolerance_ = std::ceil(error_tolerance_ * convert_ratio_);
      };

  PrsmParaPtr prsm_para_ptr_;

  std::string input_file_ext_;

  std::string output_file_ext_;

  std::string residue_mod_file_;

  int n_ = 100;

  int N_ = 10000;

  int k_ = 3;

  int max_known_mods_ = 10;

  int thread_num_ = 1;

  double convert_ratio_;

  double error_tolerance_;

  int int_tolerance_;

  int getIntTolerance() {return int_tolerance_;}
};

typedef std::shared_ptr<MCMCMng> MCMCMngPtr;

}  // namespace toppic

#endif
