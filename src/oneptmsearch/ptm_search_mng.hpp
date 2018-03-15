//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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


#ifndef PROT_PTM_SEARCH_MNG_HPP_
#define PROT_PTM_SEARCH_MNG_HPP_

#include <string>

#include "prsm/prsm_para.hpp"
#include "oneptmsearch/ps_align_para.hpp"

namespace prot {

class PtmSearchMng {
 public :
  PtmSearchMng(PrsmParaPtr prsm_para_ptr, int n_report,
               double align_max_shift,
               double align_min_shift,
               int shift_num,
               int thread_num,
               const std::string &input_file_ext,
               const std::string &output_file_ext):
      prsm_para_ptr_(prsm_para_ptr),
      n_report_(n_report),
      thread_num_(thread_num),
      input_file_ext_(input_file_ext),
      output_file_ext_(output_file_ext) {
        align_para_ptr_ = std::make_shared<PsAlignPara>(shift_num, align_max_shift, align_min_shift);
      }

  PrsmParaPtr prsm_para_ptr_;

  /* parameters for ptm search */
  int n_report_ = 1;

  int thread_num_ = 1;

  std::string input_file_ext_;
  std::string output_file_ext_;

  /* parameters for compute shift low memory */
  int ptm_fast_filter_scale_ = 100;

  int n_top_diagonals_ = 20;

  double min_double_gap_ = 0.25;

  int min_diagonal_gap_ = static_cast<int>(ptm_fast_filter_scale_ * min_double_gap_);

  /* parameters for diagonal generation */
  double extend_trunc_error_tolerance_ = 0.5;
  double align_prefix_suffix_shift_thresh_ = 300;

  PsAlignParaPtr align_para_ptr_;
};

typedef std::shared_ptr<PtmSearchMng> PtmSearchMngPtr;

} /* namespace prot */

#endif /* ONE_PTM_SEARCH_MNG_HPP_ */
