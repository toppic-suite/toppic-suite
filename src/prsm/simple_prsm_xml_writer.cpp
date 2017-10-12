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

#include <string>
#include <iostream>
#include <algorithm>
#include <vector>

#include "base/xml_dom_document.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_util.hpp"
#include "base/file_util.hpp"
#include "spec/msalign_reader.hpp"
#include "prsm/simple_prsm_xml_writer.hpp"

namespace prot {

SimplePrsmXmlWriter::SimplePrsmXmlWriter(const std::string &file_name) {
  file_.open(file_name.c_str());
  file_ << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
  file_ << "<simple_prsm_list>" << std::endl;
}

void SimplePrsmXmlWriter::close() {
  file_ << "</simple_prsm_list>" << std::endl;
  file_.close();
}

void SimplePrsmXmlWriter::write(SimplePrsmStrPtr prsm_str_ptr) {
  std::vector<std::string> strs = prsm_str_ptr->getStrVec();
  for (size_t i = 0; i < strs.size(); i++) {
    file_ << strs[i] << std::endl;
  }
}

void SimplePrsmXmlWriter::write(const SimplePrsmPtrVec &simple_prsm_ptrs) {
  for (size_t i = 0; i < simple_prsm_ptrs.size(); i++) {
    write(simple_prsm_ptrs[i]);
  }
}

void SimplePrsmXmlWriter::write(SimplePrsmPtr simple_prsm_ptr) {
  std::vector<std::string> strs = simple_prsm_ptr->toStrVec();
  for (size_t i = 0; i < strs.size(); i++) {
    file_ << strs[i] << std::endl;
  }
}

} /* namespace prot */
