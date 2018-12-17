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

#include "base/logger.hpp"
#include "base/amino_acid.hpp"
#include "base/string_util.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"
#include "base/xml_dom_util.hpp"

namespace toppic {

AminoAcid::AminoAcid(xercesc::DOMElement* element) {
  name_ = xml_dom_util::getChildValue(element, "name", 0);
  one_letter_ = xml_dom_util::getChildValue(element, "one_letter", 0);
  three_letter_ = xml_dom_util::getChildValue(element, "three_letter", 0);
  composition_ = xml_dom_util::getChildValue(element, "composition", 0);
  mono_mass_ = xml_dom_util::getDoubleChildValue(element, "mono_mass", 0);
  average_mass_ = xml_dom_util::getDoubleChildValue(element, "average_mass", 0);
}

void AminoAcid::appendNameToXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent) {
  std::string element_name = AminoAcid::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  xml_doc->addElement(element, "name", name_.c_str());
  parent->appendChild(element);
}

std::string AminoAcid::getNameFromXml(xercesc::DOMElement * element) {
  std::string name = xml_dom_util::getChildValue(element, "name", 0);
  return name;
}

}  // namespace toppic
