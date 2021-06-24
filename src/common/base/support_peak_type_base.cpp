//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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
#include "common/xml/xml_dom_parser.hpp"
#include "common/xml/xml_dom_document.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "common/base/support_peak_type_data.hpp"
#include "common/base/support_peak_type_base.hpp"

namespace toppic {

SPTypePtrVec SPTypeBase::sp_type_ptr_vec_;

SPTypePtr SPTypeBase::sp_type_ptr_N_TERM_;

void SPTypeBase::initBase() {
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (!parser) {
    LOG_ERROR("Error in parsing support peak type data!");
    exit(EXIT_FAILURE);
  }
  xercesc::MemBufInputSource mem_str((const XMLByte*)sp_type_base_data.c_str(), 
                                     sp_type_base_data.length(), 
                                     "support_peak_type_data");
  XmlDOMDocument doc(parser, mem_str);
  XmlDOMElement* parent = doc.getDocumentElement();
  std::string element_name = SupportPeakType::getXmlElementName();
  int prm_peak_type_num = xml_dom_util::getChildCount(parent, element_name.c_str());
  for (int i = 0; i < prm_peak_type_num; i++) {
    XmlDOMElement* element = xml_dom_util::getChildElement(parent, element_name.c_str(), i);
    SPTypePtr sp_type_ptr = std::make_shared<SupportPeakType>(element);
    sp_type_ptr_vec_.push_back(sp_type_ptr);
    if (sp_type_ptr->getName() == getName_N_TERM()) {
      sp_type_ptr_N_TERM_ = sp_type_ptr;
    }
  }
}

SPTypePtr SPTypeBase::getSPTypePtrByName(const std::string &name) {
  for (size_t i = 0; i < sp_type_ptr_vec_.size(); i++) {
    std::string n = sp_type_ptr_vec_[i]->getName();
    if (n == name) {
      return sp_type_ptr_vec_[i];
    }
  }
  LOG_WARN("Support peak type " << name << " cannot be found!");
  return SPTypePtr(nullptr);
}

SPTypePtr SPTypeBase::getSPTypePtrById(int id) {
  for (size_t i = 0; i < sp_type_ptr_vec_.size(); i++) {
    int n = sp_type_ptr_vec_[i]->getId();
    if (n == id) {
      return sp_type_ptr_vec_[i];
    }
  }
  LOG_WARN("Support peak id " << id << " cannot be found!");
  return SPTypePtr(nullptr);
}

}  // namespace toppic
