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


#include "base/string_util.hpp"
#include "spec/peak.hpp"

namespace prot {

Peak::Peak(double position, double intensity):
    position_(position),
    intensity_(intensity) {
    }

void Peak::appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
  std::string element_name = Peak::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  std::string str = string_util::convertToString(getPosition());
  xml_doc->addElement(element, "position", str.c_str());
  str = string_util::convertToString(getIntensity());
  xml_doc->addElement(element, "intensity", str.c_str());
  parent->appendChild(element);
}

}
