//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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
#include "ms/spec/peak_tolerance.hpp"

namespace toppic {

PeakTolerance::PeakTolerance(double ppo, bool use_min_tolerance,
                             double min_tolerance):
    ppo_(ppo), 
    use_min_tolerance_(use_min_tolerance),
    min_tolerance_(min_tolerance) {}

double PeakTolerance::compRelaxErrorTole(double m1, double m2) {
  return compStrictErrorTole(m1 + m2);
}

PeakTolerance::PeakTolerance(xercesc::DOMElement* element) {
  ppo_ = xml_dom_util::getDoubleChildValue(element, "ppo", 0);
  use_min_tolerance_ = xml_dom_util::getDoubleChildValue(element, "use_min_tolerance", 0);
  min_tolerance_ = xml_dom_util::getDoubleChildValue(element, "min_tolerance", 0);
}

double PeakTolerance::compStrictErrorTole(double mass) {
  double tolerance = mass * ppo_;
  if (use_min_tolerance_ && tolerance < min_tolerance_) {
    tolerance = min_tolerance_;
  }
  return tolerance;
}

void PeakTolerance::appendXml(XmlDOMDocument* xml_doc, 
                              xercesc::DOMElement* parent) {
  std::string element_name = PeakTolerance::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  std::string str = str_util::toString(ppo_);
  xml_doc->addElement(element, "ppo", str.c_str());
  str = str_util::toString(use_min_tolerance_);
  xml_doc->addElement(element, "use_min_tolerance", str.c_str());
  str = str_util::toString(min_tolerance_);
  xml_doc->addElement(element, "min_tolerance", str.c_str());
  parent->appendChild(element);
}

} /* namespace toppic */
