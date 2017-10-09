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


#include "base/logger.hpp"
#include "spec/msalign_util.hpp"

namespace prot {

int MsAlignUtil::countSpNum(const std::string &spectrum_file_name, SpParaPtr sp_para_ptr) {
  MsAlignReader reader(spectrum_file_name, 1, nullptr, sp_para_ptr->getSkipList());
  int cnt = 0;
  DeconvMsPtr deconv_ms_ptr;
  while ((deconv_ms_ptr = reader.getNextMs()) != nullptr) {
    cnt++;
  }
  reader.close();
  return cnt;
}

void MsAlignUtil::geneSpIndex(const std::string &spectrum_file_name, SpParaPtr sp_para_ptr) {
  int sp_num = countSpNum(spectrum_file_name, sp_para_ptr); 
  std::ofstream index_output;
  std::string index_file_name = spectrum_file_name + "_index";
  index_output.open(index_file_name.c_str(), std::ios::out);
  index_output << sp_num << std::endl;
  index_output.close();
}

int MsAlignUtil::getSpNum(const std::string &spectrum_file_name) {
  std::ifstream index_input;
  std::string index_file_name = spectrum_file_name + "_index";
  index_input.open(index_file_name.c_str(), std::ios::in);
  std::string line;
  std::getline(index_input, line);
  int sp_num = std::stoi(line);
  LOG_DEBUG("Get sp number " << sp_num);
  return sp_num;
}

}
