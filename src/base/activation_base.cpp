//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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

#include "base/activation_base.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

ActivationPtrVec ActivationBase::activation_ptr_vec_;

void ActivationBase::initBase(const std::string &file_name) {
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    std::string element_name = Activation::getXmlElementName();
    int activation_num = XmlDomUtil::getChildCount(parent, element_name.c_str());
    for (int i = 0; i < activation_num; i++) {
      xercesc::DOMElement* element
          = XmlDomUtil::getChildElement(parent, element_name.c_str(), i);
      ActivationPtr ptr(new Activation(element));
      activation_ptr_vec_.push_back(ptr);
    }
  }
}

ActivationPtr ActivationBase::getActivationPtrByName(const std::string &name) {
  for (size_t i = 0; i < activation_ptr_vec_.size(); i++) {
    std::string n = activation_ptr_vec_[i]->getName();
    if (n == name) {
      return activation_ptr_vec_[i];
    }
  }
  return ActivationPtr(nullptr);
}

ActivationPtr ActivationBase::getActivationPtrFromXml(xercesc::DOMElement * element) {
  std::string name = Activation::getNameFromXml(element);
  ActivationPtr activation_ptr = getActivationPtrByName(name);
  return activation_ptr;
}

}  // namespace prot
