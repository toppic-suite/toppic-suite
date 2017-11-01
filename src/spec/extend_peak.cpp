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
#include "spec/extend_peak.hpp"

namespace prot {

ExtendPeak::ExtendPeak(DeconvPeakPtr base_peak_ptr, 
                       double mono_mass,double score):
    Peak(mono_mass,1.0),
    base_peak_ptr_(base_peak_ptr),
    mono_mass_(mono_mass),
    score_(score),
    orig_tolerance_(0.0),
    reverse_tolerance_(0.0) {
    }

void ExtendPeak::appendXml(XmlDOMDocument* xml_doc,
                           xercesc::DOMElement* parent) {
  std::string element_name = ExtendPeak::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  std::string str = StringUtil::convertToString(getPosition());
  xml_doc->addElement(element, "position", str.c_str());
  str = StringUtil::convertToString(getIntensity());
  xml_doc->addElement(element, "intensity", str.c_str());
  str = StringUtil::convertToString(mono_mass_);
  xml_doc->addElement(element, "mono_mass", str.c_str());
  str = StringUtil::convertToString(score_);
  xml_doc->addElement(element, "score", str.c_str());
  str = StringUtil::convertToString(orig_tolerance_);
  xml_doc->addElement(element, "orig_tolerance", str.c_str());
  str = StringUtil::convertToString(reverse_tolerance_);
  xml_doc->addElement(element, "reverse_tolerance", str.c_str());
  base_peak_ptr_->appendXml(xml_doc,element);
  parent->appendChild(element);
}

} /* namespace prot */
