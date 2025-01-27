//Copyright (c) 2014 - 2025, The Trustees of Indiana University, Tulane University.
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

#ifndef TOPPIC_SEARCH_VAR_PTM_SEARCH_VAR_PTM_SEARCH_MNG_HPP_
#define TOPPIC_SEARCH_VAR_PTM_SEARCH_VAR_PTM_SEARCH_MNG_HPP_

#include <string>

#include "para/prsm_para.hpp"

namespace toppic {

class VarPtmSearchMng {
 public :
  VarPtmSearchMng(PrsmParaPtr prsm_para_ptr, int n_report,
                  const std::string &var_ptm_file_name,
                  int var_ptm_num,
                  int thread_num,
                  const std::string &input_file_ext,
                  const std::string &output_file_ext);

  void computeShifts();

  PrsmParaPtr prsm_para_ptr_;

  // parameters for variable ptm search 
  int n_report_ = 1;

  int var_ptm_num_ = 3;

  int thread_num_ = 1;

  std::string input_file_ext_;
  std::string output_file_ext_;

  // variable PTM mass shift list 
  std::vector<double> single_shift_list_;
  std::vector<int> round_single_shift_list_;
  double round_scale_ = 10000;

  // amino acid list that can be modified by the shift
  ModPtrVec2D mod_ptr_vec_2d_;

  ResiduePtrVec2D res_ptr_vec_2d_;

  // shift list for single and multiple (up to var_ptm_num) variable PTMs  
  std::vector<double> shift_list_;
  std::vector<int> round_shift_list_;

  std::vector<std::vector<int>> diag_prev_idxes_;

  std::vector<std::vector<int>> diag_prev_shift_idxes_;

  std::vector<std::vector<int>> diag_matrix_shift_idxes_;
};

typedef std::shared_ptr<VarPtmSearchMng> VarPtmSearchMngPtr;

} /* namespace toppic */

#endif /* ONE_PTM_SEARCH_MNG_HPP_ */
