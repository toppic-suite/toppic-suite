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

#include "common/util/str_util.hpp"
#include "common/xml/xml_dom_document.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "common/base/ptm.hpp"

namespace toppic {

Ptm::Ptm(const std::string &name, const std::string &abbr_name,
         double mono_mass, int unimod_id):
    name_(name),
    abbr_name_(abbr_name),
    mono_mass_(mono_mass),
    unimod_id_(unimod_id) {}

Ptm::Ptm(XmlDOMElement* element) {
  name_ = xml_dom_util::getChildValue(element, "name", 0);
  abbr_name_ = xml_dom_util::getChildValue(element, "abbreviation", 0);
  mono_mass_ = xml_dom_util::getDoubleChildValue(element, "mono_mass", 0);
  unimod_id_ = xml_dom_util::getIntChildValue(element, "unimod", 0);
}

void Ptm::appendAbbrNameToXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent) {
  std::string element_name = Ptm::getXmlElementName();
  XmlDOMElement* element = xml_doc->createElement(element_name.c_str());
  xml_doc->addElement(element, "abbreviation", abbr_name_.c_str());
  std::string str = str_util::toString(unimod_id_);
  xml_doc->addElement(element, "unimod", str.c_str());
  parent->appendChild(element);
}


void Ptm::appendAbbrNameMassToXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent) {
  std::string element_name = Ptm::getXmlElementName();
  XmlDOMElement* element = xml_doc->createElement(element_name.c_str());
  xml_doc->addElement(element, "abbreviation", abbr_name_.c_str());
  std::string str = str_util::toString(unimod_id_);
  xml_doc->addElement(element, "unimod", str.c_str());
  str = str_util::toString(mono_mass_);
  xml_doc->addElement(element, "mono_mass", str.c_str());
  parent->appendChild(element);
}

std::string Ptm::getAbbrNameFromXml(XmlDOMElement * element) {
  std::string abbr_name = xml_dom_util::getChildValue(element, "abbreviation", 0);
  return abbr_name;
}

}  // namespace toppic

