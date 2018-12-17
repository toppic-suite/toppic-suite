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


#include <string>

#include "base/support_peak_type_base.hpp"
#include "base/xml_dom_util.hpp"

namespace toppic {

SPTypePtrVec SPTypeBase::sp_type_ptr_vec_;

SPTypePtr SPTypeBase::sp_type_ptr_N_TERM_;

void SPTypeBase::initBase(const std::string &file_name) {
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    std::string element_name = SupportPeakType::getXmlElementName();
    int prm_peak_type_num = xml_dom_util::getChildCount(parent, element_name.c_str());
    for (int i = 0; i < prm_peak_type_num; i++) {
      xercesc::DOMElement* element = xml_dom_util::getChildElement(parent, element_name.c_str(), i);
      SPTypePtr sp_type_ptr = std::make_shared<SupportPeakType>(element);
      sp_type_ptr_vec_.push_back(sp_type_ptr);
      if (sp_type_ptr->getName() == getName_N_TERM()) {
        sp_type_ptr_N_TERM_ = sp_type_ptr;
      }
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
  return SPTypePtr(nullptr);
}

SPTypePtr SPTypeBase::getSPTypePtrById(int id) {
  for (size_t i = 0; i < sp_type_ptr_vec_.size(); i++) {
    int n = sp_type_ptr_vec_[i]->getId();
    if (n == id) {
      return sp_type_ptr_vec_[i];
    }
  }
  return SPTypePtr(nullptr);
}

}  // namespace toppic
