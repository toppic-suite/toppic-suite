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
#include "common/base/activation_base.hpp"

namespace toppic {

std::unordered_map<std::string, ActivationPtr> ActivationBase::activation_name_map_;

// initialize activation database 
void ActivationBase::initBase(const std::string &base_dir) {
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (!parser) {
    LOG_ERROR("Error in parsing activation data!");
    throw std::runtime_error("Error in parsing activation data!");
  }
  std::string activation_base_file_name = base_dir 
      + file_util::getFileSeparator() + "activation_base.xml";
  std::string activation_base_data = file_util::readFile(activation_base_file_name);
  xercesc::MemBufInputSource mem_str((const XMLByte*)activation_base_data.c_str(), 
                                     activation_base_data.length(), 
                                     "activation_data");
  XmlDOMDocument doc(parser, mem_str);
  XmlDOMElement* parent = doc.getDocumentElement();
  std::string element_name = Activation::getXmlElementName();
  int activation_num = xml_dom_util::getChildCount(parent, element_name.c_str());
  for (int i = 0; i < activation_num; i++) {
    XmlDOMElement* element
        = xml_dom_util::getChildElement(parent, element_name.c_str(), i);
    ActivationPtr ptr = std::make_shared<Activation>(element);
    activation_name_map_[ptr->getName()] = ptr;
  }
}

ActivationPtr ActivationBase::getActivationPtrByName(const std::string &name) {
  auto it = activation_name_map_.find(name);
  if (it == activation_name_map_.end()) {
    LOG_WARN("Activation type " << name << " is not found!");
    return nullptr;
  }
  return it->second;
}

ActivationPtr ActivationBase::getActivationPtrFromXml(XmlDOMElement* element) {
  std::string name = Activation::getNameFromXml(element);
  ActivationPtr activation_ptr = getActivationPtrByName(name);
  return activation_ptr;
}

}  // namespace toppic
