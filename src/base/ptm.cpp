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


#include <iostream>
#include <vector>
#include <memory>
#include <algorithm>
#include <string>

#include "base/ptm.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"
#include "base/xml_dom_util.hpp"
#include "base/string_util.hpp"

namespace prot {

Ptm::Ptm(xercesc::DOMElement* element) {
  name_ = XmlDomUtil::getChildValue(element, "name", 0);
  abbr_name_ = XmlDomUtil::getChildValue(element, "abbreviation", 0);
  mono_mass_ = XmlDomUtil::getDoubleChildValue(element, "mono_mass", 0);
  unimod_id_ = XmlDomUtil::getIntChildValue(element, "unimod", 0);
}

void Ptm::appendAbbrNameToXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent) {
  std::string element_name = Ptm::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  xml_doc->addElement(element, "abbreviation", abbr_name_.c_str());
  xml_doc->addElement(element, "unimod", StringUtil::convertToString(unimod_id_).c_str());
  parent->appendChild(element);
}

std::string Ptm::getAbbrNameFromXml(xercesc::DOMElement * element) {
  std::string abbr_name = XmlDomUtil::getChildValue(element, "abbreviation", 0);
  return abbr_name;
}

}  // namespace prot

