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

#include "common/base/mod_util.hpp"
#include "filter/mng/var_ptm_filter_mng.hpp"
#include "search/varptmsearch/var_ptm_search_mng.hpp"

namespace toppic {

VarPtmSearchMng::VarPtmSearchMng(PrsmParaPtr prsm_para_ptr, int n_report,
                                 const std::string &var_ptm_file_name,
                                 int var_ptm_num,
                                 int thread_num,
                                 const std::string &input_file_ext,
                                 const std::string &output_file_ext):
  prsm_para_ptr_(prsm_para_ptr),
  n_report_(n_report),
  var_ptm_num_(var_ptm_num), 
  thread_num_(thread_num),
  input_file_ext_(input_file_ext),
  output_file_ext_(output_file_ext) {
    single_shift_list_ = mod_util::readModTxtToShiftList(var_ptm_file_name);
    shift_list_ = VarPtmFilterMng::computeShifts(single_shift_list_, var_ptm_num_);
  }

} /* namespace toppic */
