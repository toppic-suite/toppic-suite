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


#ifndef PROT_PRSM_PRSM_STR_COMBINE_HPP_
#define PROT_PRSM_PRSM_STR_COMBINE_HPP_

#include <vector>
#include <string>
#include <map>

#include "base/proteoform.hpp"
#include "base/fasta_reader.hpp"
#include "prsm/prsm.hpp"
#include "prsm/prsm_para.hpp"
#include "prsm/prsm_xml_writer.hpp"

namespace prot {

class PrsmStrCombine {
 public:
  PrsmStrCombine(const std::string &spec_file_name, 
                 const std::vector<std::string> &in_file_exts,
                 const std::string &out_file_ext,
                 int top_num);

  PrsmStrCombine(const std::string &spec_file_name, 
                 const std::string &in_file_ext,
                 int in_num,
                 const std::string &out_file_ext, 
                 int top_num);

  void process(bool norm);

  void process() {process(false);}

 private:
  std::string spec_file_name_;
  std::vector<std::string> input_file_exts_;
  std::string output_file_ext_;
  int top_num_;
};

typedef std::shared_ptr<PrsmStrCombine> PrsmStrCombinePtr;
} /* namespace prot */

#endif 
