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

#include "base/string_util.hpp"
#include "base/xml_dom_util.hpp"
#include "spec/env_peak.hpp"

namespace prot {

EnvPeak::EnvPeak(double mz, double intensity):
    Peak(mz, intensity) {
      idx_ = EnvPeak::getNonExistPeakIdx();
    }

EnvPeak::EnvPeak(xercesc::DOMElement* element):
    Peak(xml_dom_util::getDoubleChildValue(element, "position", 0),
         xml_dom_util::getDoubleChildValue(element, "intensity", 0)) {
      idx_ = xml_dom_util::getIntChildValue(element, "index", 0);
    }

void EnvPeak::appendXml(XmlDOMDocument* xml_doc, 
                        xercesc::DOMElement* parent) {
  std::string element_name = EnvPeak::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  std::string str = string_util::convertToString(getPosition());
  xml_doc->addElement(element, "position", str.c_str());
  str = string_util::convertToString(getIntensity());
  xml_doc->addElement(element, "intensity", str.c_str());
  str = string_util::convertToString(idx_);
  xml_doc->addElement(element, "index", str.c_str());
  parent->appendChild(element);
}

bool EnvPeak::isExist() {
  if (idx_ != getNonExistPeakIdx()) {
    return true;
  } else {
    return false;
  }
}

}  // namespace prot

