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

#include "feature/single_charge_feature.hpp"

namespace toppic {

SingleChargeFeature::SingleChargeFeature(int charge,
                                         double time_begin, double time_end,
                                         int scan_begin, int scan_end,
                                         double intensity, int env_num): 
    charge_(charge),
    time_begin_(time_begin),
    time_end_(time_end),
    scan_begin_(scan_begin),
    scan_end_(scan_end),
    intensity_(intensity),
    env_num_(env_num) {
    }

SingleChargeFeature::SingleChargeFeature(XmlDOMElement* element) {
  charge_ = xml_dom_util::getIntChildValue(element, "charge", 0);
  time_begin_ = xml_dom_util::getDoubleChildValue(element, "time_begin", 0);
  time_end_ = xml_dom_util::getDoubleChildValue(element, "time_end", 0);
  scan_begin_ = xml_dom_util::getIntChildValue(element, "scan_begin", 0);
  scan_end_ = xml_dom_util::getIntChildValue(element, "scan_end", 0);
  intensity_ = xml_dom_util::getDoubleChildValue(element, "intensity", 0);
  env_num_ = xml_dom_util::getIntChildValue(element, "envelope_num", 0);
}

void SingleChargeFeature::appendToXml(XmlDOMDocument* xml_doc, XmlDOMElement *parent) {
  std::string element_name = SingleChargeFeature::getXmlElementName();
  XmlDOMElement* element = xml_doc->createElement(element_name.c_str());
  std::string str = str_util::toString(charge_);
  xml_doc->addElement(element, "charge", str.c_str());
  str = str_util::toString(time_begin_);
  xml_doc->addElement(element, "time_begin", str.c_str());
  str = str_util::toString(time_end_);
  xml_doc->addElement(element, "time_end", str.c_str());
  str = str_util::toString(scan_begin_);
  xml_doc->addElement(element, "scan_begin", str.c_str());
  str = str_util::toString(scan_end_);
  xml_doc->addElement(element, "scan_end", str.c_str());
  str = str_util::toString(intensity_);
  xml_doc->addElement(element, "intensity", str.c_str());
  str = str_util::toString(env_num_);
  xml_doc->addElement(element, "envelope_num", str.c_str());
  parent->appendChild(element);
}

}
