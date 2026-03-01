//Copyright (c) 2014 - 2025, The Trustees of Indiana University, Tulane University.
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

#include <stdexcept>
#include <string>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "common/xml/xml_dom_parser.hpp"
#include "common/xml/xml_dom_document.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "common/base/neutral_loss_base.hpp"

namespace toppic {

NeutralLossPtrVec NeutralLossBase::neutral_loss_ptr_vec_;
std::unordered_map<std::string, NeutralLossPtr> NeutralLossBase::neutral_loss_name_map_;

NeutralLossPtr NeutralLossBase::neutral_loss_ptr_NONE_;

void NeutralLossBase::initBase(const std::string &base_dir) {
  XmlDOMParser * parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (!parser) {
    LOG_ERROR("Error in parsing neutral loss data!");
    throw std::runtime_error("Error in parsing neutral loss data!");
  }

  std::string neutral_loss_base_file_name = base_dir 
      + file_util::getFileSeparator() + "neutral_loss_base.xml";
  std::string neutral_loss_base_data = file_util::readFile(neutral_loss_base_file_name);
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
    neutral_loss_name_map_[neutral_loss_ptr->getName()] = neutral_loss_ptr;
    if (neutral_loss_ptr->getName() == getName_NONE()) {
      neutral_loss_ptr_NONE_ = neutral_loss_ptr;
    }
  }
}

NeutralLossPtr NeutralLossBase::getNeutralLossPtrByName(const std::string &name) {
  auto it = neutral_loss_name_map_.find(name);
  if (it == neutral_loss_name_map_.end()) {
    LOG_WARN("Neutral loss " << name << " cannot be found!");
    return nullptr;
  }
  return it->second;
}

} /* namespace toppic */
