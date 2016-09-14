// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#ifndef PROT_PTM_SEARCH_MNG_HPP_
#define PROT_PTM_SEARCH_MNG_HPP_

#include "prsm/prsm_para.hpp"
#include "oneptmsearch/ps_align_para.hpp"

namespace prot {

class PtmSearchMng {
 public :
  PtmSearchMng(PrsmParaPtr prsm_para_ptr, int n_report, 
               double align_max_shift, int shift_num, 
               int thread_num,
               const std::string &input_file_ext, 
               const std::string &output_file_ext) {
    prsm_para_ptr_ = prsm_para_ptr;
    n_report_ = n_report;
    thread_num_ = thread_num;
    input_file_ext_ = input_file_ext;
    output_file_ext_ = output_file_ext;
    align_para_ptr_ = PsAlignParaPtr(new PsAlignPara(shift_num, align_max_shift));
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
  double min_double_gap=0.25;
  int min_diagonal_gap_ = (int)(ptm_fast_filter_scale_ * min_double_gap);

  /* parameters for diagonal generation */
  double extend_trunc_error_tolerance_ = 0.5;
  double align_prefix_suffix_shift_thresh_ = 300;

  PsAlignParaPtr align_para_ptr_;
};

typedef std::shared_ptr<PtmSearchMng> PtmSearchMngPtr;

} /* namespace prot */

#endif /* ONE_PTM_SEARCH_MNG_HPP_ */
