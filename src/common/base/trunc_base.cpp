//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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

#include <string>

#include "common/util/logger.hpp"
#include "common/xml/xml_dom_document.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "common/base/trunc_data.hpp"
#include "common/base/trunc_base.hpp"

namespace toppic {

TruncPtrVec TruncBase::trunc_ptr_vec_;

void TruncBase::initBase() {
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (!parser) {
    LOG_ERROR("Error in parsing truncation data!");
    exit(EXIT_FAILURE);
  }
  xercesc::MemBufInputSource mem_str((const XMLByte*)trunc_base_data.c_str(), 
                                     trunc_base_data.length(), 
                                     "truncation_data");
  XmlDOMDocument doc(parser, mem_str);
  XmlDOMElement* parent = doc.getDocumentElement();
  std::string element_name = Trunc::getXmlElementName();
  int trunc_num = xml_dom_util::getChildCount(parent, element_name.c_str());
  for (int i = 0; i < trunc_num; i++) {
    XmlDOMElement* element = xml_dom_util::getChildElement(parent, element_name.c_str(), i);
    TruncPtr trunc_ptr = std::make_shared<Trunc>(element);
    trunc_ptr_vec_.push_back(trunc_ptr);
  }
}

TruncPtr TruncBase::getTruncPtrByName(const std::string &name) {
  for (size_t i = 0; i < trunc_ptr_vec_.size(); i++) {
    std::string n = trunc_ptr_vec_[i]->getName();
    if (n == name) {
      return trunc_ptr_vec_[i];
    }
  }
  LOG_ERROR("Truncation " << name << " cannot be found!");
  return TruncPtr(nullptr);
}

TruncPtr TruncBase::getTruncPtrFromXml(XmlDOMElement * element) {
  std::string name = Trunc::getNameFromXml(element);
  TruncPtr trunc_ptr = getTruncPtrByName(name);
  return trunc_ptr;
}

}  // namespace toppic

