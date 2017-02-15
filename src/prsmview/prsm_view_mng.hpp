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


#ifndef PROT_PRSM_VIEW_MNG_HPP_
#define PROT_PRSM_VIEW_MNG_HPP_

#include "base/file_util.hpp"
#include "prsm/prsm_para.hpp"

namespace prot {

class PrsmViewMng {
 public:
  PrsmViewMng(PrsmParaPtr prsm_para_ptr, std::string exec_dir);

  PrsmParaPtr prsm_para_ptr_;

  std::string html_path_;
  std::string xml_path_;
  std::string executive_dir_;

  int decimal_point_num_ = 2;
  int precise_point_num_ = 4;
  double min_mass_;
  size_t cnt_;
  size_t num_files_;
};

typedef std::shared_ptr<PrsmViewMng> PrsmViewMngPtr;

inline PrsmViewMng::PrsmViewMng(PrsmParaPtr prsm_para_ptr, std::string exec_dir) {
  prsm_para_ptr_ = prsm_para_ptr;
  std::string spectrum_file_name = prsm_para_ptr_->getSpectrumFileName();

  xml_path_ = FileUtil::basename(spectrum_file_name) + "_xml";
  html_path_ = FileUtil::basename(spectrum_file_name) + "_html";
  executive_dir_ = exec_dir;
  min_mass_ = prsm_para_ptr_->getSpParaPtr()->getMinMass();
  cnt_ = 0;
}

} /* namespace prot */


#endif /* PROT_PRSM_VIEW_MNG_HPP_ */
