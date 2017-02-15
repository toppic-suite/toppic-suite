// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#ifndef PROT_ZERO_PTM_SEARCH_MNG_HPP_
#define PROT_ZERO_PTM_SEARCH_MNG_HPP_

#include <string>

#include "prsm/prsm_para.hpp"

namespace prot {

class ZeroPtmSearchMng {
 public:
  ZeroPtmSearchMng(PrsmParaPtr prsm_para_ptr, 
                   const std::string & input_file_ext,
                   const std::string & output_file_ext): 
      prsm_para_ptr_(prsm_para_ptr), 
      input_file_ext_(input_file_ext),
      output_file_ext_(output_file_ext) {
      }

  PrsmParaPtr prsm_para_ptr_;

  /** parameters for zero ptm search */

  /** zero ptm fast filtering */
  int zero_ptm_filter_result_num_ = 20;
  /** number of reported Prsms for each spectrum */
  int report_num_ = 1;

  /** recalibration is used in ZeroPtmSlowMatch */
  bool   do_recalibration_ = false;
  double recal_ppo_ = 0.000015; // 15 ppm
  bool   ms_one_ms_two_same_recal_ = true;

  std::string input_file_ext_;

  std::string output_file_ext_;
};

typedef std::shared_ptr<ZeroPtmSearchMng> ZeroPtmSearchMngPtr;

} /* namespace_prot */

#endif
