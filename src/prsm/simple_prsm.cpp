// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
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


#include <iostream>
#include <cmath>
#include <algorithm>

#include "base/logger.hpp"
#include "base/proteoform_factory.hpp"
#include "base/proteoform.hpp"
#include "base/string_util.hpp"
#include "base/xml_dom_document.hpp"
#include "base/xml_dom_util.hpp"

#include "prsm/simple_prsm.hpp"

namespace prot {

SimplePrsm::SimplePrsm(MsHeaderPtr header_ptr, int spectrum_num, 
                       ProteoformPtr proteo_ptr,int score): 
  spectrum_num_(spectrum_num),
  //proteo_ptr_(proteo_ptr),
  score_(score) {
  spectrum_id_ = header_ptr->getId();
  spectrum_scan_ = header_ptr->getScansString();
  precursor_id_ = header_ptr->getPrecId();
  prec_mass_ = header_ptr->getPrecMonoMass();
  seq_name_ = proteo_ptr->getSeqName();
  seq_desc_ = proteo_ptr->getSeqDesc();
  prot_mass_ = proteo_ptr->getResSeqPtr()->getSeqMass();
}

SimplePrsm::SimplePrsm(xercesc::DOMElement* element){
  spectrum_id_ = XmlDomUtil::getIntChildValue(element, "spectrum_id", 0);
  spectrum_scan_ = XmlDomUtil::getChildValue(element, "spectrum_scan", 0);
  precursor_id_ = XmlDomUtil::getIntChildValue(element, "precursor_id", 0);
  prec_mass_ = XmlDomUtil::getDoubleChildValue(element, "precursor_mass", 0);
  spectrum_num_ = XmlDomUtil::getDoubleChildValue(element, "spectrum_number", 0);
  seq_name_ = XmlDomUtil::getChildValue(element, "sequence_name", 0);
  seq_desc_ = XmlDomUtil::getChildValue(element, "sequence_desc", 0);
  prot_mass_ = XmlDomUtil::getDoubleChildValue(element, "proteoform_mass", 0);
  score_ = XmlDomUtil::getDoubleChildValue(element, "score", 0);
  // n trunc shifts
  xercesc::DOMElement* n_shift_list_element= XmlDomUtil::getChildElement(element, "n_trunc_shift_list",0);
  int n_shift_num = XmlDomUtil::getChildCount(n_shift_list_element, "shift");
  //LOG_DEBUG("n shift _num " << n_shift_num);
  for(int i=0; i<n_shift_num; i++) {
    double shift = XmlDomUtil::getDoubleChildValue(n_shift_list_element, "shift", i);
    n_trunc_shifts_.push_back(shift);
  }
  //c trunc shifts
  xercesc::DOMElement* c_shift_list_element= XmlDomUtil::getChildElement(element, "c_trunc_shift_list",0);
  int c_shift_num = XmlDomUtil::getChildCount(c_shift_list_element, "shift");
  //LOG_DEBUG("c shift _num " << c_shift_num);
  for(int i=0; i<c_shift_num; i++) {
    double shift = XmlDomUtil::getDoubleChildValue(c_shift_list_element, "shift", i);
    c_trunc_shifts_.push_back(shift);
  }
}

xercesc::DOMElement* SimplePrsm::toXml(XmlDOMDocument* xml_doc){
  std::string element_name = SimplePrsm::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  std::string str = StringUtil::convertToString(spectrum_id_);
  xml_doc->addElement(element, "spectrum_id", str.c_str());
  xml_doc->addElement(element, "spectrum_scan", spectrum_scan_.c_str());
  str = StringUtil::convertToString(precursor_id_);
  xml_doc->addElement(element, "precursor_id", str.c_str());
  str = StringUtil::convertToString(prec_mass_);
  xml_doc->addElement(element, "precursor_mass", str.c_str());
  str = StringUtil::convertToString(spectrum_num_);
  xml_doc->addElement(element, "spectrum_number", str.c_str());
  xml_doc->addElement(element, "sequence_name", seq_name_.c_str());
  xml_doc->addElement(element, "sequence_desc", seq_desc_.c_str());
  str = StringUtil::convertToString(prot_mass_);
  xml_doc->addElement(element, "proteoform_mass", str.c_str());
  str = StringUtil::convertToString(score_);
  xml_doc->addElement(element, "score", str.c_str());

  xercesc::DOMElement* n_shift_list = xml_doc->createElement("n_trunc_shift_list");
  for(size_t i=0; i<n_trunc_shifts_.size(); i++) {
    str = StringUtil::convertToString(n_trunc_shifts_[i]);
    xml_doc->addElement(n_shift_list, "shift", str.c_str());
  }
  element->appendChild(n_shift_list);

  xercesc::DOMElement* c_shift_list = xml_doc->createElement("c_trunc_shift_list");
  for(size_t i=0; i<c_trunc_shifts_.size(); i++) {
    str = StringUtil::convertToString(c_trunc_shifts_[i]);
    xml_doc->addElement(c_shift_list, "shift", str.c_str());
  }
  element->appendChild(c_shift_list);
  return element;
}

void SimplePrsm::setCTruncShifts(const std::vector<double> &c_term_shifts) {
  for (size_t i = 0; i < c_term_shifts.size(); i++) {
    double shift = prec_mass_ - (prot_mass_ + c_term_shifts[i]);
    c_trunc_shifts_.push_back(shift);
  }
}

} /* namespace prot */
