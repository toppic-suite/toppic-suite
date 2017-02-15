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


#include <limits>

#include "base/logger.hpp"
#include "prsm/prsm_util.hpp"
#include "prsm/prsm_str.hpp"

namespace prot {

PrsmStr::PrsmStr(const std::vector<std::string> &str_vec) {
  str_vec_ = str_vec;
  std::string line = PrsmUtil::getXmlLine(str_vec_, "<spectrum_id>");
  spectrum_id_ = std::stoi(PrsmUtil::getValueStr(line));
  line = PrsmUtil::getXmlLine(str_vec_, "<seq_name>");
  seq_name_ = PrsmUtil::getValueStr(line);
  line = PrsmUtil::getXmlLine(str_vec_, "<match_fragment_num>");
  match_frag_num_ = std::stod(PrsmUtil::getValueStr(line));
  line = PrsmUtil::getXmlLine(str_vec_, "<norm_match_fragment_num>");
  norm_match_frag_num_ = std::stod(PrsmUtil::getValueStr(line));
  line = PrsmUtil::getXmlLine(str_vec_, "<e_value>"); 
  if (line == "") { 
    e_value_ = 0.0;
  } else {
    e_value_ = std::stod(PrsmUtil::getValueStr(line));
  }
  line = PrsmUtil::getXmlLine(str_vec_, "<fdr>"); 
  fdr_ = std::stod(PrsmUtil::getValueStr(line));
  line = PrsmUtil::getXmlLine(str_vec_, "<proteoform_fdr>"); 
  proteoform_fdr_ = std::stod(PrsmUtil::getValueStr(line));
  line = PrsmUtil::getXmlLine(str_vec_, "<species_id>");
  species_id_ = std::stoi(PrsmUtil::getValueStr(line));
  //LOG_DEBUG("spectrum id " << spectrum_id_ << " match num " << match_frag_num_);
}

int getXmlLineIndex(const std::vector<std::string> &str_vec,
                    const std::string &property) {
  for (size_t i = 0; i < str_vec.size(); i++) {
    size_t found = str_vec[i].find(property);
    if (found != std::string::npos) {
      return i;
    }
  }
  return -1;
}

void PrsmStr::setFdr(double fdr) {
  int i = getXmlLineIndex(str_vec_, "fdr");
  str_vec_[i] = "<fdr>" + std::to_string(fdr) + "</fdr>";
}

void PrsmStr::setProteoformFdr(double proteoform_fdr) {
  int i = getXmlLineIndex(str_vec_, "proteoform_fdr");
  str_vec_[i] = "<proteoform_fdr>" + std::to_string(proteoform_fdr) + "</proteoform_fdr>";
}

void PrsmStr::setId(int id) {
  int i = getXmlLineIndex(str_vec_, "prsm_id");
  str_vec_[i] = "<prsm_id>" + std::to_string(id) + "</prsm_id>";
}

}

