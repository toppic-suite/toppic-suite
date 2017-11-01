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


#ifndef PROT_POSSION_MNG_HPP_
#define PROT_POSSION_MNG_HPP_

#include "prsm/prsm_para.hpp"

namespace prot {

class PoissonMng {
 public:
  PoissonMng(PrsmParaPtr prsm_para_ptr, int shift_num, double max_ptm_mass,
             std::string input_file_ext, std::string output_file_ext) {
    prsm_para_ptr_ = prsm_para_ptr,
    max_ptm_mass_ = max_ptm_mass;
    unexpected_shift_num_ = shift_num;
    input_file_ext_ = input_file_ext;
    output_file_ext_ = output_file_ext;
  }

  PrsmParaPtr prsm_para_ptr_;

  /** parameters  */
  double max_ptm_mass_ = 1000000;

  int comp_evalue_min_match_frag_num_ = 4;

  // number of mass shift
  int unexpected_shift_num_ = 2;
  double convert_ratio_ = 274.335215;
  double max_prec_mass_ = 51000;

  int shift_penalty = 4;

  std::string input_file_ext_;
  std::string output_file_ext_;
};

typedef std::shared_ptr<PoissonMng> PoissonMngPtr;

}
#endif
