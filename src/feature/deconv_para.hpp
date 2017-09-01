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


#ifndef PROT_FEATURE_DECONV_PARA_HPP_
#define PROT_FEATURE_DECONV_PARA_HPP_

#include <string>
#include <memory>
#include <map>
#include <vector>

namespace prot {

enum InputType {INPUT_MGF, INPUT_MZXML};

class DeconvPara {
 public:
  DeconvPara(std::map<std::string, std::string> &arguments);

  void setDataFileName(const std::string & file_name) {data_file_name_ = file_name;}

  std::string getDataFileName() {return data_file_name_;}

  int setInputType (const std::string & format);

  std::string data_file_name_;
  std::string exec_dir_;
  InputType input_type_;

  bool refine_prec_mass_;
  bool missing_level_one_;

  int max_charge_;
  double max_mass_;
  double tolerance_;
  double sn_ratio_;
  bool keep_unused_peaks_;

  bool output_multiple_mass_ = false; 

  double prec_window_;
};

typedef std::shared_ptr<DeconvPara> DeconvParaPtr;

}
#endif
