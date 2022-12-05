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
                                 bool use_approx_spec,
                                 const std::string &output_file_ext): 
  prsm_para_ptr_(prsm_para_ptr),
  index_file_para_(index_file_para),
  var_ptm_num_(var_ptm_num),
  thread_num_(thread_num),
  use_approx_spec_(use_approx_spec),
  output_file_ext_(output_file_ext) {
    single_shift_list_ = mod_util::readModTxtToShiftList(var_ptm_file_name);
    // round single shift list
    for (size_t i = 0; i < single_shift_list_.size(); i++) {
      int round_shift = std::lround(single_shift_list_[i] * round_scale_);
      round_single_shift_list_.push_back(round_shift);
    }
    LOG_DEBUG("Number of single shifts:" << round_single_shift_list_.size());
    // use integer to avoid errors in double addition
    round_shift_list_ = computeShifts(round_single_shift_list_, var_ptm_num_);
    for (size_t i = 0; i < round_shift_list_.size(); i++) {
      double shift = round_shift_list_[i] /round_scale_; 
      shift_list_.push_back(shift);
    }

    // convert to integer again using filter scaling factor
    LOG_DEBUG("Number of shifts:" << round_shift_list_.size());
    for (size_t i = 0; i < shift_list_.size(); i++) {
      LOG_DEBUG("Shifts:" << i << " " << shift_list_[i]);
      int int_shift = std::round(shift_list_[i] * filter_scale_); 
      int_shift_list_.push_back(int_shift);
    }
  }

std::vector<int> VarPtmFilterMng::computeShifts(std::vector<int> &round_single_shift_list,
                                                int var_ptm_num) {
  std::vector<int> round_shift_list;
  round_shift_list.push_back(0);
  for (int i = 0; i < var_ptm_num; i++) {
    size_t cur_list_len = round_shift_list.size();
    for (size_t j = 0; j < cur_list_len; j++) {
      for (size_t k = 0; k < round_single_shift_list.size(); k++) {
        int new_shift = round_shift_list[j] + round_single_shift_list[k];
        // if the new shift is not in the list
        if (std::find(round_shift_list.begin(), round_shift_list.end(), new_shift)
            == round_shift_list.end()) {
          round_shift_list.push_back(new_shift);
        }
      }
    }
  }
  return round_shift_list;
}
}
