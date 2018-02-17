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

#include "base/mass_shift_type.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

const MassShiftTypePtr MassShiftType::INPUT            = std::make_shared<MassShiftType>(1, "Input");

const MassShiftTypePtr MassShiftType::FIXED            = std::make_shared<MassShiftType>(2, "Fixed");

const MassShiftTypePtr MassShiftType::PROTEIN_VARIABLE = std::make_shared<MassShiftType>(3, "Protein variable");

const MassShiftTypePtr MassShiftType::VARIABLE         = std::make_shared<MassShiftType>(4, "Variable");

const MassShiftTypePtr MassShiftType::UNEXPECTED       = std::make_shared<MassShiftType>(5, "Unexpected");

void MassShiftType::appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent) {
  std::string element_name = MassShiftType::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  xml_doc->addElement(element, "name", name_.c_str());
  parent->appendChild(element);
}

MassShiftTypePtr MassShiftType::getChangeTypePtrFromXml(xercesc::DOMElement * element) {
  std::string name = xml_dom_util::getChildValue(element, "name", 0);
  if (name == MassShiftType::INPUT->getName()) {
    return MassShiftType::INPUT;
  }
  if (name == MassShiftType::FIXED->getName()) {
    return MassShiftType::FIXED;
  }
  if (name == MassShiftType::PROTEIN_VARIABLE->getName()) {
    return MassShiftType::PROTEIN_VARIABLE;
  }
  if (name == MassShiftType::VARIABLE->getName()) {
    return MassShiftType::VARIABLE;
  }
  if (name == MassShiftType::UNEXPECTED->getName()) {
    return MassShiftType::UNEXPECTED;
  }
  return nullptr;
}

}  // namespace prot
