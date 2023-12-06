//Copyright (c) 2014 - 2023, The Trustees of Indiana University.
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

#include "ms/msmap/ms_map_peak.hpp"

namespace toppic {

MsMapPeak::MsMapPeak(PeakPtr peak):
  Peak(peak->getPosition(), peak->getIntensity()) {
    ori_inte_ = peak->getIntensity();
    neighbor_ = false;
  }

std::string MsMapPeak::getString() {
  return "Pos: " + std::to_string(getPosition()) + " " +
    "Inte: " + std::to_string(getIntensity()) + " " +
    "Neighbor: " + std::to_string(neighbor_) + "\n";
}

void MsMapPeak::appendToXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent) {
  std::string element_name = "ms_map_peak";
  XmlDOMElement* element = xml_doc->createElement(element_name.c_str());
  std::string str = str_util::toString(getPosition());
  xml_doc->addElement(element, "position", str.c_str());
  str = str_util::toString(getIntensity());
  xml_doc->addElement(element, "intensity", str.c_str());
  str = str_util::toString(ori_inte_);
  xml_doc->addElement(element, "ori_intensity", str.c_str());
  parent->appendChild(element);
}


}
