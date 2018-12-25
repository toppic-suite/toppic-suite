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

#include "common/xml/xml_dom_parser.hpp"
#include "common/xml/xml_dom_document.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "common/base/neutral_loss_base.hpp"

namespace toppic {

NeutralLossPtrVec NeutralLossBase::neutral_loss_ptr_vec_;

NeutralLossPtr NeutralLossBase::neutral_loss_ptr_NONE_;

void NeutralLossBase::initBase(const std::string &file_name) {
  XmlDOMParser * parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    toppic::XmlDOMDocument doc(parser, file_name.c_str());
    XmlDOMElement* parent = doc.getDocumentElement();
    std::string element_name = NeutralLoss::getXmlElementName();
    int neutral_loss_num = xml_dom_util::getChildCount(parent, element_name.c_str());
    for (int i = 0; i < neutral_loss_num; i++) {
      XmlDOMElement* element = xml_dom_util::getChildElement(parent, element_name.c_str(), i);
      NeutralLossPtr neutral_loss_ptr = std::make_shared<NeutralLoss>(element);
      neutral_loss_ptr_vec_.push_back(neutral_loss_ptr);
      if (neutral_loss_ptr->getName() == getName_NONE()) {
        neutral_loss_ptr_NONE_ = neutral_loss_ptr;
      }
    }
  }
}

NeutralLossPtr NeutralLossBase::getNeutralLossPtrByName(const std::string &name) {
  for (size_t i = 0; i < neutral_loss_ptr_vec_.size(); i++) {
    std::string n = neutral_loss_ptr_vec_[i]->getName();
    if (n == name) {
      return neutral_loss_ptr_vec_[i];
    }
  }
  return NeutralLossPtr(nullptr);
}

} /* namespace toppic */
