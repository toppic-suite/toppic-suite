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

#include "common/util/str_util.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "common/xml/xml_dom_document.hpp"
#include "ms/spec/env_peak.hpp"

namespace toppic {

EnvPeak::EnvPeak(double mz, double intensity):
    Peak(mz, intensity) {
      idx_ = EnvPeak::getNonExistPeakIdx();
    }

EnvPeak::EnvPeak(double mz, double intensity, int idx):
      Peak(mz, intensity),
      idx_(idx) {}

EnvPeak::EnvPeak(EnvPeakPtr peak_ptr):
      Peak(peak_ptr->getPosition(), peak_ptr->getIntensity()),
      idx_(peak_ptr->getIdx()) {}


bool EnvPeak::cmpPosInc(const EnvPeakPtr &a, const EnvPeakPtr &b) {
  return a->getPosition() < b->getPosition();
}

bool EnvPeak::cmpInteInc(const EnvPeakPtr &a, const EnvPeakPtr &b) {
  return a->getIntensity() < b->getIntensity();
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
  std::string str = str_util::toString(getPosition());
  xml_doc->addElement(element, "position", str.c_str());
  str = str_util::toString(getIntensity());
  xml_doc->addElement(element, "intensity", str.c_str());
  str = str_util::toString(idx_);
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

}  // namespace toppic

