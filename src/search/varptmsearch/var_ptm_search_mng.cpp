// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "search/varptmsearch/var_ptm_search_mng.hpp"

#include <cmath>
#include <cstdlib>

#include "common/base/mod_util.hpp"
#include "common/base/ptm.hpp"
#include "common/base/residue.hpp"
#include "common/util/logger.hpp"

namespace toppic {

VarPtmSearchMng::VarPtmSearchMng(const PrsmParaPtr& prsm_para_ptr, int n_report,
                                 const std::string& var_ptm_file_name,
                                 int var_ptm_num, int thread_num,
                                 const std::string& input_file_ext,
                                 const std::string& output_file_ext)
    : prsm_para_ptr_(prsm_para_ptr),
      n_report_(n_report),
      var_ptm_num_(var_ptm_num),
      thread_num_(thread_num),
      input_file_ext_(input_file_ext),
      output_file_ext_(output_file_ext) {
  ModPtrVec mod_ptr_list = mod_util::readAnywhereModTxt(var_ptm_file_name);
  for (size_t i = 0; i < mod_ptr_list.size(); i++) {
    double shift = mod_ptr_list[i]->getShift();
    int pos = -1;
    for (size_t j = 0; j < single_shift_list_.size(); j++) {
      if (single_shift_list_[j] == shift) {
        pos = j;
        break;
      }
    }
    // if the shift is not in the list
    if (pos == -1) {
      single_shift_list_.push_back(shift);
      ModPtrVec cur_mod_ptr_vec;
      cur_mod_ptr_vec.push_back(mod_ptr_list[i]);
      mod_ptr_vec_2d_.push_back(cur_mod_ptr_vec);

      ResiduePtrVec cur_res_ptr_vec;
      cur_res_ptr_vec.push_back(mod_ptr_list[i]->getOriResiduePtr());
      res_ptr_vec_2d_.push_back(cur_res_ptr_vec);
    } else {
      mod_ptr_vec_2d_[pos].push_back(mod_ptr_list[i]);
      res_ptr_vec_2d_[pos].push_back(mod_ptr_list[i]->getOriResiduePtr());
    }
    // the modifications of the same PTM in one list
    const PtmPtr& ptm_ptr = mod_ptr_list[i]->getModResiduePtr()->getPtmPtr();
    int ptm_pos = -1;
    for (size_t j = 0; j < ptm_mod_ptr_vec_2d_.size(); j++) {
      if (ptm_mod_ptr_vec_2d_[j][0]->getModResiduePtr()->getPtmPtr()->isSame(
              ptm_ptr)) {
        ptm_pos = j;
        break;
      }
    }
    if (ptm_pos == -1) {
      ptm_mod_ptr_vec_2d_.push_back(ModPtrVec{mod_ptr_list[i]});
    } else {
      ptm_mod_ptr_vec_2d_[ptm_pos].push_back(mod_ptr_list[i]);
    }
  }
  n_term_mod_ptr_vec_ = mod_util::readNTermModTxt(var_ptm_file_name);
  // get round single shifts
  // round single shift list
  for (size_t i = 0; i < single_shift_list_.size(); i++) {
    int round_shift = std::lround(single_shift_list_[i] * round_scale_);
    round_single_shift_list_.push_back(round_shift);
  }

  computeShifts();
}

const ModPtrVec& VarPtmSearchMng::getPtmModPtrVec(
    const ModPtr& mod_ptr) const {
  const PtmPtr& ptm_ptr = mod_ptr->getModResiduePtr()->getPtmPtr();
  for (size_t i = 0; i < ptm_mod_ptr_vec_2d_.size(); i++) {
    if (ptm_mod_ptr_vec_2d_[i][0]->getModResiduePtr()->getPtmPtr()->isSame(
            ptm_ptr)) {
      return ptm_mod_ptr_vec_2d_[i];
    }
  }
  LOG_ERROR("The PTM " << ptm_ptr->getAbbrName()
                       << " is not a variable PTM!");
  exit(EXIT_FAILURE);
}

void VarPtmSearchMng::computeShifts() {
  round_shift_list_.push_back(0);
  // empty previous indexes for shift 0
  std::vector<int> empty_idxes;
  diag_prev_idxes_.push_back(empty_idxes);
  std::vector<int> empty_shift_idxes;
  diag_prev_shift_idxes_.push_back(empty_shift_idxes);

  for (int i = 0; i < var_ptm_num_; i++) {
    size_t cur_list_len = round_shift_list_.size();
    for (size_t j = 0; j < cur_list_len; j++) {
      for (size_t k = 0; k < round_single_shift_list_.size(); k++) {
        int new_shift = round_shift_list_[j] + round_single_shift_list_[k];
        LOG_DEBUG("new shift " << new_shift);
        int pos = -1;
        for (size_t p = 0; p < round_shift_list_.size(); p++) {
          if (round_shift_list_[p] == new_shift) {
            pos = p;
            break;
          }
        }
        // if the new shift is not in the list
        if (pos == -1) {
          round_shift_list_.push_back(new_shift);
          std::vector<int> cur_idxes;
          cur_idxes.push_back(j);
          diag_prev_idxes_.push_back(cur_idxes);
          std::vector<int> shift_idxes;
          shift_idxes.push_back(k);
          diag_prev_shift_idxes_.push_back(shift_idxes);
        } else {
          diag_prev_idxes_[pos].push_back(j);
          diag_prev_shift_idxes_[pos].push_back(k);
        }
      }
    }
  }

  // get double shifts
  LOG_DEBUG("Number of shifts:" << round_shift_list_.size());
  for (size_t i = 0; i < round_shift_list_.size(); i++) {
    double shift = round_shift_list_[i] / round_scale_;
    shift_list_.push_back(shift);
    LOG_DEBUG("Shifts:" << i << " " << shift_list_[i]);
  }

  // generate diag_matrix_shift_idxes
  for (size_t i = 0; i < diag_prev_idxes_.size(); i++) {
    std::vector<int> matrix_row(diag_prev_idxes_.size(), -1);
    for (size_t j = 0; j < diag_prev_idxes_[i].size(); j++) {
      int prev = diag_prev_idxes_[i][j];
      int shift_idx = diag_prev_shift_idxes_[i][j];
      matrix_row[prev] = shift_idx;
    }
    diag_matrix_shift_idxes_.push_back(matrix_row);
  }

  /*
  for (size_t i = 0; i < diag_prev_idxes_.size(); i++) {
    for (size_t j = 0; j < diag_prev_idxes_.size(); j++) {
      std::cout << diag_matrix_shift_idxes_[i][j] << " ";
    }
    std::cout << std::endl;
  }
  */
}

} /* namespace toppic */
