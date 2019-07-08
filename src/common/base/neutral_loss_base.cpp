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

#include "common/util/logger.hpp"
#include "common/xml/xml_dom_parser.hpp"
#include "common/xml/xml_dom_document.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "common/base/neutral_loss_data.hpp"
#include "common/base/neutral_loss_base.hpp"

namespace toppic {

NeutralLossPtrVec NeutralLossBase::neutral_loss_ptr_vec_;

NeutralLossPtr NeutralLossBase::neutral_loss_ptr_NONE_;

void NeutralLossBase::initBase() {
  XmlDOMParser * parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (!parser) {
    LOG_ERROR("Error in parsing neutral loss data!");
    exit(EXIT_FAILURE);
  }

  xercesc::MemBufInputSource mem_str((const XMLByte*)neutral_loss_base_data.c_str(), 
                                     neutral_loss_base_data.length(), 
                                     "neutral_loss_data");
  XmlDOMDocument doc(parser, mem_str); 
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
