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

#include "util/str_util.hpp"
#include "xml/xml_dom_util.hpp"
#include "spec/deconv_peak.hpp"

namespace toppic {

DeconvPeak::DeconvPeak(xercesc::DOMElement* element):
    Peak(xml_dom_util::getDoubleChildValue(element, "position", 0),
         xml_dom_util::getDoubleChildValue(element, "intensity", 0)) {
      id_ = xml_dom_util::getIntChildValue(element, "id", 0);
      charge_ = xml_dom_util::getIntChildValue(element, "charge", 0);
      score_ = xml_dom_util::getDoubleChildValue(element, "score", 0);
    }

void DeconvPeak::appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent) {
  std::string element_name = DeconvPeak::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  std::string str = str_util::toString(getPosition());
  xml_doc->addElement(element, "position", str.c_str());
  str = str_util::toString(getIntensity());
  xml_doc->addElement(element, "intensity", str.c_str());
  str = str_util::toString(id_);
  xml_doc->addElement(element, "id", str.c_str());
  str = str_util::toString(charge_);
  xml_doc->addElement(element, "charge", str.c_str());
  str = str_util::toString(score_);
  xml_doc->addElement(element, "score", str.c_str());
  parent->appendChild(element);
}

}  // namespace toppic

