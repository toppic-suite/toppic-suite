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

#include "common/xml/xml_dom_util.hpp"
#include "common/xml/xml_dom_document.hpp"
#include "seq/alter_type.hpp"

namespace toppic {

const AlterTypePtr AlterType::INPUT 
    = std::make_shared<AlterType>(1, "Input");

const AlterTypePtr AlterType::FIXED 
    = std::make_shared<AlterType>(2, "Fixed");

const AlterTypePtr AlterType::PROTEIN_VARIABLE 
    = std::make_shared<AlterType>(3, "Protein variable");

const AlterTypePtr AlterType::VARIABLE 
    = std::make_shared<AlterType>(4, "Variable");

const AlterTypePtr AlterType::UNEXPECTED
    = std::make_shared<AlterType>(5, "Unexpected");

AlterType::AlterType(int id, std::string name): 
    id_(id), 
    name_(name) {}

void AlterType::appendXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent) {
  std::string element_name = AlterType::getXmlElementName();
  XmlDOMElement* element = xml_doc->createElement(element_name.c_str());
  xml_doc->addElement(element, "name", name_.c_str());
  parent->appendChild(element);
}

AlterTypePtr AlterType::getTypePtrFromXml(XmlDOMElement * element) {
  std::string name = xml_dom_util::getChildValue(element, "name", 0);
  if (name == AlterType::INPUT->getName()) {
    return AlterType::INPUT;
  }
  if (name == AlterType::FIXED->getName()) {
    return AlterType::FIXED;
  }
  if (name == AlterType::PROTEIN_VARIABLE->getName()) {
    return AlterType::PROTEIN_VARIABLE;
  }
  if (name == AlterType::VARIABLE->getName()) {
    return AlterType::VARIABLE;
  }
  if (name == AlterType::UNEXPECTED->getName()) {
    return AlterType::UNEXPECTED;
  }
  return nullptr;
}

}  // namespace toppic
