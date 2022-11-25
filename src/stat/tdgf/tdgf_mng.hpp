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

#ifndef TOPPIC_STAT_TDGF_TDGF_MNG_HPP_
#define TOPPIC_STAT_TDGF_TDGF_MNG_HPP_

#include <string>

#include "para/prsm_para.hpp"

namespace toppic {

class TdgfMng {
 public:
  TdgfMng(PrsmParaPtr prsm_para_ptr, int shift_num, double max_ptm_mass, bool use_gf, 
          int var_ptm_type_num, int thread_num, const std::string &input_file_ext, 
          const std::string & output_file_ext);

  std::string input_file_ext_;

  std::string output_file_ext_;

  PrsmParaPtr prsm_para_ptr_;

  // Prsm filter 
  double comp_evalue_min_match_frag_num_ = 4.0;

  bool use_gf_ = true;

  int var_ptm_type_num_ = 0;

  // Max ptm mass is used in the function for counting sequence numbers
  double max_ptm_mass_ = 1000000;

  // Do tdgf computation if poisson report evalue > 10^-8 
  // or match frag num < 25 
  double computation_evalue_cutoff = 0.00000001;

  int computation_frag_num_cutoff = 25;

  // dp table 
  // number of mass shift
  int unexpected_shift_num_ = 2;

  int thread_num_ = 1;

  double convert_ratio_ = 27.4335215;

  double max_prec_mass_ = 100000;  

  int max_table_height_ = 128;

  int min_height_ = 10;
};

typedef std::shared_ptr<TdgfMng> TdgfMngPtr;

}  // namespace toppic

#endif
