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


#include <xercesc/framework/MemBufInputSource.hpp>

#include "base/logger.hpp"
#include "base/string_util.hpp"
#include "prsm/simple_prsm_reader.hpp"

namespace prot {

SimplePrsmReader::SimplePrsmReader(const std::string &file_name) {
  input_.open(file_name.c_str(), std::ios::in);
}

std::vector<std::string> SimplePrsmReader::readOnePrsmLines() {
  std::string line;
  std::vector<std::string> line_list;
  while (std::getline(input_, line)) {
    line = StringUtil::trim(line);
    if (line ==  "<simple_prsm>") {
      line_list.push_back(line);
    }
    else if (line == "</simple_prsm>") {
      if (line_list.size() != 0) {
        line_list.push_back(line);
      }
      return line_list;
    }
    else if (line == "") {
      continue;
    }
    else {
      if (line_list.size() > 0) {
        line_list.push_back(line);
      }
    }
  }
  return line_list;
}

SimplePrsmStrPtr SimplePrsmReader::readOnePrsmStr() {
  std::vector<std::string> prsm_str_vec = readOnePrsmLines();
  //LOG_DEBUG("prsm str size " << prsm_str_vec.size());
  if (prsm_str_vec.size() == 0) {
    return SimplePrsmStrPtr(nullptr);
  }
  return SimplePrsmStrPtr(new SimplePrsmStr(prsm_str_vec));
}

SimplePrsmPtr SimplePrsmReader::readOnePrsm() {
  std::vector<std::string> prsm_str_vec = readOnePrsmLines();
  if (prsm_str_vec.size() == 0) {
    return SimplePrsmPtr(nullptr);
  }
  std::string prsm_str = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>";
  for (size_t i = 0; i < prsm_str_vec.size(); i++) {
    prsm_str += prsm_str_vec[i];
  }
  xercesc::MemBufInputSource prsm_buf(
      (const XMLByte*)prsm_str.c_str(), prsm_str.size(), "prsm_str (in memory)");

  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  SimplePrsmPtr ptr;
  if(parser){
    XmlDOMDocument doc(parser, prsm_buf);
    xercesc::DOMElement* root = doc.getDocumentElement();
    ptr = SimplePrsmPtr(new SimplePrsm(root));
  }
  //LOG_DEBUG("simple prsm spectrum id " << ptr->getSpectrumId() << " seq name " << ptr->getSeqName());
  return ptr;
}


void SimplePrsmReader::close() {
  input_.close();
}

SimplePrsmPtrVec SimplePrsmReader::readSimplePrsms(const std::string &file_name){
  SimplePrsmPtrVec result_ptrs;
  SimplePrsmReader reader(file_name);
  SimplePrsmPtr prsm_ptr = reader.readOnePrsm();
  while (prsm_ptr != nullptr) {
    result_ptrs.push_back(prsm_ptr);
    prsm_ptr = reader.readOnePrsm();
  }
  return result_ptrs;
}


} /* namespace prot */
