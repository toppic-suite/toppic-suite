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


#include "base/file_util.hpp"
#include "prsm/prsm_str_combine.hpp"
#include "prsm/prsm_reader.hpp"
#include "prsm/prsm_str.hpp"

namespace prot {

PrsmStrCombine::PrsmStrCombine(const std::string &spec_file_name, 
                               const std::vector<std::string> &in_file_exts,
                               const std::string &out_file_ext, 
                               int top_num) {
  input_file_exts_ = in_file_exts;
  output_file_ext_ = out_file_ext;
  spec_file_name_ = spec_file_name;
  top_num_ = top_num;
}

PrsmStrCombine::PrsmStrCombine(const std::string &spec_file_name, 
                               const std::string &in_file_ext,
                               int in_num,
                               const std::string &out_file_ext, 
                               int top_num) {
  output_file_ext_ = out_file_ext;
  spec_file_name_ = spec_file_name;
  top_num_ = top_num;
  for (int i = 0; i < in_num; i ++) {
    std::string ext = in_file_ext + "_" + std::to_string(i);
    input_file_exts_.push_back(ext);
  }
}

void PrsmStrCombine::process(bool norm) {
  size_t input_num = input_file_exts_.size();
  std::string base_name = FileUtil::basename(spec_file_name_); 
  // open files
  PrsmReaderPtrVec reader_ptrs;
  PrsmStrPtrVec prsm_str_ptrs;
  for (size_t i = 0; i < input_num; i++) {
    std::string input_file_name = base_name + "." + input_file_exts_[i]; 
    PrsmReaderPtr reader_ptr(new PrsmReader(input_file_name));
    PrsmStrPtr str_ptr = reader_ptr->readOnePrsmStr();
    reader_ptrs.push_back(reader_ptr);
    prsm_str_ptrs.push_back(str_ptr);
  }
  PrsmXmlWriter writer(base_name +"."+output_file_ext_);
  
  // combine
  int spec_id = 0;
  bool finish = false;
  while (!finish) {
    finish = true;
    PrsmStrPtrVec cur_str_ptrs;
    for (size_t i = 0; i < input_num; i++) {
      if (prsm_str_ptrs[i] != nullptr) {
        finish = false;
        while (prsm_str_ptrs[i] != nullptr && prsm_str_ptrs[i]->getSpectrumId() == spec_id) {
          cur_str_ptrs.push_back(prsm_str_ptrs[i]);
          prsm_str_ptrs[i] = reader_ptrs[i]->readOnePrsmStr();
        }
      }
    }
    if (cur_str_ptrs.size() > 0) {
      if (!norm) {
        std::sort(cur_str_ptrs.begin(),cur_str_ptrs.end(), PrsmStr::cmpMatchFragmentDec);
      }
      else {
        std::sort(cur_str_ptrs.begin(),cur_str_ptrs.end(), PrsmStr::cmpNormMatchFragmentDec);
      }
      for (int i = 0; i < top_num_; i++) {
        if (i >= (int)cur_str_ptrs.size()) {
          break;
        }
        writer.write(cur_str_ptrs[i]);
      }
    }
    spec_id++;
  }
  
  // close files
  for (size_t i = 0; i < input_num; i++) {
    reader_ptrs[i]->close();
  }

  writer.close();
}

} /* namespace prot */
