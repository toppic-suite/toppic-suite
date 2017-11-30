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

#include "base/change_type.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

const ChangeTypePtr ChangeType::INPUT = std::make_shared<ChangeType>(1, "Input");

const ChangeTypePtr ChangeType::FIXED = std::make_shared<ChangeType>(2, "Fixed");

const ChangeTypePtr ChangeType::PROTEIN_VARIABLE = std::make_shared<ChangeType>(3, "Protein variable");

const ChangeTypePtr ChangeType::VARIABLE = std::make_shared<ChangeType>(4, "Variable");

const ChangeTypePtr ChangeType::UNEXPECTED = std::make_shared<ChangeType>(5, "Unexpected");

void ChangeType::appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent) {
  std::string element_name = ChangeType::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  xml_doc->addElement(element, "name", name_.c_str());
  parent->appendChild(element);
}

ChangeTypePtr ChangeType::getChangeTypePtrFromXml(xercesc::DOMElement * element) {
  std::string name = xml_dom_util::getChildValue(element, "name", 0);
  if (name == ChangeType::INPUT->getName()) {
    return ChangeType::INPUT;
  }
  if (name == ChangeType::FIXED->getName()) {
    return ChangeType::FIXED;
  }
  if (name == ChangeType::PROTEIN_VARIABLE->getName()) {
    return ChangeType::PROTEIN_VARIABLE;
  }
  if (name == ChangeType::VARIABLE->getName()) {
    return ChangeType::VARIABLE;
  }
  if (name == ChangeType::UNEXPECTED->getName()) {
    return ChangeType::UNEXPECTED;
  }
  return nullptr;
}

}  // namespace prot
