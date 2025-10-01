// Copyright (c) 2014 - 2025, The Trustees of Indiana University, Tulane University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
#include <algorithm>

#include "ms/feature/single_charge_feature.hpp"

#include "common/util/str_util.hpp"
#include "common/xml/xml_dom_util.hpp"

namespace toppic {

SingleChargeFeature::SingleChargeFeature(
		int charge, double time_begin, double time_end, int scan_begin,
		int scan_end, double intensity, int env_num, int spec_id_begin,
		int spec_id_end, double mass, std::vector<double> xic_inte,
		std::vector<double> envelope_mass, 
		std::vector<double> aggregate_envelope_inte,
		double mono_mz, double average_mz,
		double refer_mz,
		std::vector<int> scan_list,
		std::vector<double> rt_list,
		std::vector<double> intensity_sum_list,
		std::vector<double> max_intensity_list)
    : charge_(charge),
      time_begin_(time_begin),
      time_end_(time_end),
      scan_begin_(scan_begin),
      scan_end_(scan_end),
      intensity_(intensity),
      env_num_(env_num),
      spec_id_begin_(spec_id_begin),
      spec_id_end_(spec_id_end),
      mass_(mass),
      xic_inte_(xic_inte),
      envelope_mass_(envelope_mass),
      aggregate_envelope_inte_(aggregate_envelope_inte),
      mono_mz_(mono_mz),
      average_mz_(average_mz),
      refer_mz_(refer_mz),
      scan_list_(scan_list),
      rt_list_(rt_list),
      intensity_sum_list_(intensity_sum_list),
      max_intensity_list_(max_intensity_list) {
} 

// the constructor needs to be UPDATED. 
SingleChargeFeature::SingleChargeFeature(XmlDOMElement* element) {
  charge_ = xml_dom_util::getIntChildValue(element, "charge", 0);
  time_begin_ = xml_dom_util::getDoubleChildValue(element, "time_begin", 0);
  time_end_ = xml_dom_util::getDoubleChildValue(element, "time_end", 0);
  scan_begin_ = xml_dom_util::getIntChildValue(element, "scan_begin", 0);
  scan_end_ = xml_dom_util::getIntChildValue(element, "scan_end", 0);
  intensity_ = xml_dom_util::getDoubleChildValue(element, "intensity", 0);
  env_num_ = xml_dom_util::getIntChildValue(element, "envelope_num", 0);
  mono_mz_ = xml_dom_util::getDoubleChildValue(element, "mono_mz", 0);
  average_mz_ = xml_dom_util::getDoubleChildValue(element, "average_mz", 0);
  refer_mz_ = xml_dom_util::getDoubleChildValue(element, "refer_mz", 0);
  scan_list_.clear();
  rt_list_.clear();
  intensity_sum_list_.clear();
  max_intensity_list_.clear();
  XmlDOMElement* env_list_element = xml_dom_util::getChildElement(element, "envelope_list", 0);
  if (env_list_element) {
    for (int i = 0; i < env_num_; i++) {
      XmlDOMElement* env_element = xml_dom_util::getChildElement(env_list_element, "envelope", i);
      if (env_element) {
        scan_list_.push_back(xml_dom_util::getIntChildValue(env_element, "scan", 0));
        rt_list_.push_back(xml_dom_util::getDoubleChildValue(env_element, "rt", 0));
        intensity_sum_list_.push_back(xml_dom_util::getDoubleChildValue(env_element, "intensity_sum", 0));
        max_intensity_list_.push_back(xml_dom_util::getDoubleChildValue(env_element, "max_intensity", 0));
      }
    }
  }
}

void SingleChargeFeature::appendToXml(XmlDOMDocument* xml_doc,
                                      XmlDOMElement* parent) {
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
  str = str_util::toString(mono_mz_);
  xml_doc->addElement(element, "mono_mz", str.c_str());
  str = str_util::toString(average_mz_);
  xml_doc->addElement(element, "average_mz", str.c_str());
  str = str_util::toString(refer_mz_);
  xml_doc->addElement(element, "refer_mz", str.c_str());
  double max_inte = *std::max_element(max_intensity_list_.begin(), max_intensity_list_.end());
  str = str_util::toString(max_inte);
  xml_doc->addElement(element, "max_intensity", str.c_str());
  XmlDOMElement* env_list_element = xml_doc->createElement("envelope_list"); 
  element->appendChild(env_list_element);
  for (int i = 0; i < env_num_; i++) {
    XmlDOMElement* env_element = xml_doc->createElement("envelope");
    str = str_util::toString(scan_list_[i]);
    xml_doc->addElement(env_element, "scan", str.c_str());
    str = str_util::fixedToString(rt_list_[i], 2);
    xml_doc->addElement(env_element, "rt", str.c_str());
    str = str_util::fixedToString(intensity_sum_list_[i], 2);
    xml_doc->addElement(env_element, "intensity_sum", str.c_str());
    str = str_util::fixedToString(max_intensity_list_[i], 2);
    xml_doc->addElement(env_element, "max_intensity", str.c_str());
    env_list_element->appendChild(env_element);  
  }
  parent->appendChild(element);
}

}  // namespace toppic
