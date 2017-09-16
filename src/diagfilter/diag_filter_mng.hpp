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


#ifndef PROT_DIAG_FILTER_MNG_HPP_
#define PROT_DIAG_FILTER_MNG_HPP_

#include <string>

#include "prsm/prsm_para.hpp"

namespace prot {

class DiagFilterMng {
 public:
  DiagFilterMng(PrsmParaPtr prsm_para_ptr,
                int filtering_result_num,
                int thread_num,
                const std::string &output_file_ext,
                const std::string & residueModFileName = "",
                int var_num = 0):
      prsm_para_ptr_(prsm_para_ptr),
      filter_result_num_(filtering_result_num),
      thread_num_(thread_num),
      output_file_ext_(output_file_ext),
      residueModFileName_(residueModFileName),
      var_num_(var_num) {}

  PrsmParaPtr prsm_para_ptr_;

  /** parameters for fast filteration */
  int max_proteoform_mass_ = 20000;

  // Candidate protein number for each spectrum
  size_t filter_result_num_ = 20;
  int db_block_size_ = 5000000;
  int filter_scale_ = 100;
  int thread_num_ = 1;

  std::string output_file_ext_;
  std::string residueModFileName_;

  int var_num_;
};

typedef std::shared_ptr<DiagFilterMng> DiagFilterMngPtr;

}  // namespace prot

#endif /* DIAG_FILTER_MNG_HPP_ */
