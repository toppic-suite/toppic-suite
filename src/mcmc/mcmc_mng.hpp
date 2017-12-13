//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#ifndef PROT_MCMC_MNG_HPP_
#define PROT_MCMC_MNG_HPP_

#include "prsm/prsm_para.hpp"

namespace prot {

class MCMCMng {
 public:
  MCMCMng(PrsmParaPtr prsm_para_ptr, 
          const std::string & input_file_ext, 
          const std::string & output_file_ext,
          const std::string & residue_mod_file):
      prsm_para_ptr_(prsm_para_ptr),
      input_file_ext_(input_file_ext),
      output_file_ext_(output_file_ext),
      residue_mod_file_(residue_mod_file) {};

  PrsmParaPtr prsm_para_ptr_;
  std::string input_file_ext_;
  std::string output_file_ext_;
  std::string residue_mod_file_;

  int n_ = 300;

  int N_ = 1000000;

  int k_ = 30;

  double mass_limit_ = 150.0;
};

typedef std::shared_ptr<MCMCMng> MCMCMngPtr;

}  // namespace prot

#endif
