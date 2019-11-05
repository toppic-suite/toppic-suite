//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#include <xercesc/framework/MemBufInputSource.hpp>

#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "prsm/simple_prsm_reader.hpp"

namespace toppic {

SimplePrsmReader::SimplePrsmReader(const std::string &file_name) {
  input_.open(file_name.c_str(), std::ios::in);
}

std::vector<std::string> SimplePrsmReader::readOnePrsmLines() {
  std::string line;
  std::vector<std::string> line_list;
  while (std::getline(input_, line)) {
    str_util::trim(line);
    if (line ==  "<simple_prsm>") {
      line_list.push_back(line);
    } else if (line == "</simple_prsm>") {
      if (line_list.size() != 0) {
        line_list.push_back(line);
      }
      return line_list;
    } else if (line == "") {
      continue;
    } else {
      if (line_list.size() > 0) {
        line_list.push_back(line);
      }
    }
  }
  return line_list;
}

SimplePrsmStrPtr SimplePrsmReader::readOnePrsmStr() {
  std::vector<std::string> prsm_str_vec = readOnePrsmLines();
  // LOG_DEBUG("prsm str size " << prsm_str_vec.size());
  if (prsm_str_vec.size() == 0) {
    return SimplePrsmStrPtr(nullptr);
  }
  return std::make_shared<SimplePrsmStr>(prsm_str_vec);
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
  if (parser) {
    XmlDOMDocument doc(parser, prsm_buf);
    XmlDOMElement* root = doc.getDocumentElement();
    ptr = std::make_shared<SimplePrsm>(root);
  }
  // LOG_DEBUG("simple prsm spectrum id " << ptr->getSpectrumId() << " seq name " << ptr->getSeqName());
  return ptr;
}


void SimplePrsmReader::close() {
  input_.close();
}

SimplePrsmPtrVec SimplePrsmReader::readSimplePrsms(const std::string &file_name) {
  SimplePrsmPtrVec result_ptrs;
  SimplePrsmReader reader(file_name);
  SimplePrsmPtr prsm_ptr = reader.readOnePrsm();
  while (prsm_ptr != nullptr) {
    result_ptrs.push_back(prsm_ptr);
    prsm_ptr = reader.readOnePrsm();
  }
  return result_ptrs;
}

} /* namespace toppic */
