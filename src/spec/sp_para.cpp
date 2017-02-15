// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include "base/activation_base.hpp"
#include "base/string_util.hpp"
#include "base/xml_dom_util.hpp"
#include "spec/sp_para.hpp"

namespace prot {

SpPara::SpPara(int min_peak_num,double min_mass,
               double min_extend_mass, 
               const std::vector<double> &ext_offsets,
               PeakTolerancePtr peak_tolerance_ptr,
               ActivationPtr activation_ptr): 
    min_peak_num_(min_peak_num),
    min_mass_(min_mass),
    extend_min_mass_(min_extend_mass),
    ext_offsets_(ext_offsets),
    peak_tolerance_ptr_(peak_tolerance_ptr),
    activation_ptr_(activation_ptr) {
    }

SpPara::SpPara(xercesc::DOMElement* element){
  min_peak_num_ = XmlDomUtil::getIntChildValue(element,"min_peak_num",0);
  min_mass_ = XmlDomUtil::getDoubleChildValue(element,"min_mass",0);
  extend_min_mass_ = XmlDomUtil::getDoubleChildValue(element,"extend_min_mass",0);
  xercesc::DOMElement* list_element 
      = XmlDomUtil::getChildElement(element, "extend_offset_list", 0);
  int offset_num =  XmlDomUtil::getChildCount(list_element, "extend_offset");
  for (int i = 0; i < offset_num; i++) {
    double offset = XmlDomUtil::getDoubleChildValue(list_element, "extend_offset", i);
    ext_offsets_.push_back(offset);
  }
  std::string element_name = PeakTolerance::getXmlElementName();
  xercesc::DOMElement* pt_element = XmlDomUtil::getChildElement(element,element_name.c_str(),0);
  peak_tolerance_ptr_ = PeakTolerancePtr(new PeakTolerance(pt_element));

  element_name = Activation::getXmlElementName();
  xercesc::DOMElement* ac_element = XmlDomUtil::getChildElement(element,element_name.c_str(),0);
  activation_ptr_ = ActivationBase::getActivationPtrFromXml(ac_element);
}

void SpPara::appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent) {
  std::string element_name = SpPara::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  std::string str = StringUtil::convertToString(min_peak_num_);
  xml_doc->addElement(element, "min_peak_num", str.c_str());
  str = StringUtil::convertToString(min_mass_);
  xml_doc->addElement(element, "min_mass", str.c_str());
  str = StringUtil::convertToString(extend_min_mass_).c_str();
  xml_doc->addElement(element, "extend_min_mass", str.c_str());
  xercesc::DOMElement* list_element 
      = xml_doc->createElement("extend_offset_list");
  element->appendChild(list_element);
  for (size_t i = 0; i < ext_offsets_.size(); i++) {
    str = StringUtil::convertToString(ext_offsets_[i]);
    xml_doc->addElement(list_element, "extend_offset", str.c_str());
  }
  element->appendChild(list_element);
  peak_tolerance_ptr_->appendXml(xml_doc, element);
  activation_ptr_->appendNameToXml(xml_doc,element);
  parent->appendChild(element); 
}

} /* namespace prot */
