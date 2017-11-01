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
#include "base/xml_dom_util.hpp"
#include "spec/deconv_peak.hpp"

namespace prot {

DeconvPeak::DeconvPeak (int id, double mono_mass, 
                        double intensity, int charge):
    Peak (mono_mass, intensity),
    id_(id),
    charge_(charge),
    score_(1.0) {
    }

DeconvPeak::DeconvPeak(xercesc::DOMElement* element):
    Peak (XmlDomUtil::getDoubleChildValue(element,"position",0), 
          XmlDomUtil::getDoubleChildValue(element,"intensity",0)) {
      id_ = XmlDomUtil::getIntChildValue(element,"id",0);
      charge_ = XmlDomUtil::getIntChildValue(element,"charge",0);
      score_ = XmlDomUtil::getDoubleChildValue(element,"score",0);
    }

void DeconvPeak::appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
  std::string element_name = DeconvPeak::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  std::string str = StringUtil::convertToString(getPosition());
  xml_doc->addElement(element, "position", str.c_str());
  str = StringUtil::convertToString(getIntensity());
  xml_doc->addElement(element, "intensity", str.c_str());
  str = StringUtil::convertToString(id_);
  xml_doc->addElement(element, "id", str.c_str());
  str = StringUtil::convertToString(charge_);
  xml_doc->addElement(element, "charge", str.c_str());
  str = StringUtil::convertToString(score_);
  xml_doc->addElement(element, "score", str.c_str());
  parent->appendChild(element);
}

}

