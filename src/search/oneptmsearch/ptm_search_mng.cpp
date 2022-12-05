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

#include "search/oneptmsearch/ptm_search_mng.hpp"

namespace toppic {

PtmSearchMng::PtmSearchMng(PrsmParaPtr prsm_para_ptr, int n_report,
                           double align_min_shift,
                           double align_max_shift,
                           int shift_num,
                           int thread_num,
                           const std::string &input_file_ext,
                           const std::string &output_file_ext):
  prsm_para_ptr_(prsm_para_ptr),
  n_report_(n_report),
  thread_num_(thread_num),
  input_file_ext_(input_file_ext),
  output_file_ext_(output_file_ext) {
    align_para_ptr_ = std::make_shared<PsAlignPara>(shift_num, align_min_shift, align_max_shift);
  }

} /* namespace toppic */
