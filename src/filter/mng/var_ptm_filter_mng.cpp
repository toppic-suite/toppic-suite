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

#include "common/util/logger.hpp"
#include "common/base/mod_util.hpp"
#include "filter/mng/var_ptm_filter_mng.hpp"

namespace toppic {

VarPtmFilterMng::VarPtmFilterMng(PrsmParaPtr prsm_para_ptr,
                                 const std::string &index_file_para,
                                 const std::string &var_ptm_file_name, 
                                 int var_ptm_num,
                                 int thread_num,
                                 const std::string &output_file_ext): 
  prsm_para_ptr_(prsm_para_ptr),
  index_file_para_(index_file_para),
  var_ptm_num_(var_ptm_num),
  thread_num_(thread_num),
  output_file_ext_(output_file_ext) {
    single_shift_list_ = mod_util::readModTxtToShiftList(var_ptm_file_name);
    shift_list_.push_back(0);
    for (int i = 0; i <= var_ptm_num; i++) {
      size_t cur_list_len = shift_list_.size();
      for (size_t j = 0; j < cur_list_len; j++) {
        for (size_t k = 0; k < single_shift_list_.size(); k++) {
          double new_shift = shift_list_[j] + single_shift_list_[k];
          // if the new shift is not in the list
          if (std::find(shift_list_.begin(), shift_list_.end(), new_shift) != shift_list_.end()) {
            shift_list_.push_back(new_shift);
          }
        }
      }
    }
    LOG_DEBUG("Number of shifts:" << shift_list_.size());
    for (size_t i = 0; i < shift_list_.size(); i++) {
      LOG_DEBUG("Shifts:" << i << " " << shift_list_[i]);
    }
  }
}
