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


#include <iostream>
#include <cmath>
#include <algorithm>
#include <string>
#include <vector>

#include "base/logger.hpp"
#include "base/proteoform_factory.hpp"
#include "base/proteoform.hpp"
#include "base/string_util.hpp"
#include "base/xml_dom_document.hpp"
#include "base/xml_dom_util.hpp"

#include "prsm/simple_prsm.hpp"

namespace prot {

SimplePrsm::SimplePrsm(MsHeaderPtr header_ptr, int spectrum_num,
                       ProteoformPtr proteo_ptr, int score):
    spectrum_num_(spectrum_num),
    score_(score) {
      spectrum_id_ = header_ptr->getId();
      spectrum_scan_ = header_ptr->getScansString();
      precursor_id_ = header_ptr->getPrecId();
      prec_mass_ = header_ptr->getPrecMonoMass();
      seq_name_ = proteo_ptr->getSeqName();
      seq_desc_ = proteo_ptr->getSeqDesc();
      prot_mass_ = proteo_ptr->getResSeqPtr()->getSeqMass();
    }

SimplePrsm::SimplePrsm(xercesc::DOMElement* element) {
  spectrum_id_ = xml_dom_util::getIntChildValue(element, "spectrum_id", 0);
  spectrum_scan_ = xml_dom_util::getChildValue(element, "spectrum_scan", 0);
  precursor_id_ = xml_dom_util::getIntChildValue(element, "precursor_id", 0);
  prec_mass_ = xml_dom_util::getDoubleChildValue(element, "precursor_mass", 0);
  spectrum_num_ = xml_dom_util::getDoubleChildValue(element, "spectrum_number", 0);
  seq_name_ = xml_dom_util::getChildValue(element, "sequence_name", 0);
  seq_desc_ = xml_dom_util::getChildValue(element, "sequence_desc", 0);
  prot_mass_ = xml_dom_util::getDoubleChildValue(element, "proteoform_mass", 0);
  score_ = xml_dom_util::getDoubleChildValue(element, "score", 0);
  // n trunc shifts
  xercesc::DOMElement* n_shift_list_element
      = xml_dom_util::getChildElement(element, "n_trunc_shift_list", 0);
  int n_shift_num = xml_dom_util::getChildCount(n_shift_list_element, "shift");
  // LOG_DEBUG("n shift _num " << n_shift_num);
  for (int i = 0; i < n_shift_num; i++) {
    double shift = xml_dom_util::getDoubleChildValue(n_shift_list_element, "shift", i);
    n_trunc_shifts_.push_back(shift);
  }
  // c trunc shifts
  xercesc::DOMElement* c_shift_list_element
      = xml_dom_util::getChildElement(element, "c_trunc_shift_list", 0);
  int c_shift_num = xml_dom_util::getChildCount(c_shift_list_element, "shift");
  // LOG_DEBUG("c shift _num " << c_shift_num);
  for (int i = 0; i < c_shift_num; i++) {
    double shift = xml_dom_util::getDoubleChildValue(c_shift_list_element, "shift", i);
    c_trunc_shifts_.push_back(shift);
  }
}

xercesc::DOMElement* SimplePrsm::toXml(XmlDOMDocument* xml_doc) {
  std::string element_name = SimplePrsm::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  std::string str = string_util::convertToString(spectrum_id_);
  xml_doc->addElement(element, "spectrum_id", str.c_str());
  xml_doc->addElement(element, "spectrum_scan", spectrum_scan_.c_str());
  str = string_util::convertToString(precursor_id_);
  xml_doc->addElement(element, "precursor_id", str.c_str());
  str = string_util::convertToString(prec_mass_);
  xml_doc->addElement(element, "precursor_mass", str.c_str());
  str = string_util::convertToString(spectrum_num_);
  xml_doc->addElement(element, "spectrum_number", str.c_str());
  xml_doc->addElement(element, "sequence_name", seq_name_.c_str());
  xml_doc->addElement(element, "sequence_desc", seq_desc_.c_str());
  str = string_util::convertToString(prot_mass_);
  xml_doc->addElement(element, "proteoform_mass", str.c_str());
  str = string_util::convertToString(score_);
  xml_doc->addElement(element, "score", str.c_str());

  xercesc::DOMElement* n_shift_list = xml_doc->createElement("n_trunc_shift_list");
  for (size_t i = 0; i < n_trunc_shifts_.size(); i++) {
    str = string_util::convertToString(n_trunc_shifts_[i]);
    xml_doc->addElement(n_shift_list, "shift", str.c_str());
  }
  element->appendChild(n_shift_list);

  xercesc::DOMElement* c_shift_list = xml_doc->createElement("c_trunc_shift_list");
  for (size_t i = 0; i < c_trunc_shifts_.size(); i++) {
    str = string_util::convertToString(c_trunc_shifts_[i]);
    xml_doc->addElement(c_shift_list, "shift", str.c_str());
  }
  element->appendChild(c_shift_list);
  return element;
}

std::vector<std::string> SimplePrsm::toStrVec() {
  std::vector<std::string> str_vec;
  str_vec.push_back("<simple_prsm>");
  str_vec.push_back("<spectrum_id>" + std::to_string(getSpectrumId()) + "</spectrum_id>");
  str_vec.push_back("<spectrum_scan>" + getSpectrumScan() + "</spectrum_scan>");
  str_vec.push_back("<precursor_id>" + std::to_string(getPrecursorId()) + "</precursor_id>");
  str_vec.push_back("<precursor_mass>" + std::to_string(getPrecMass()) + "</precursor_mass>");
  str_vec.push_back("<spectrum_number>" + std::to_string(getSpectrumNum()) + "</spectrum_number>");
  str_vec.push_back("<sequence_name>" + getSeqName() + "</sequence_name>");
  str_vec.push_back("<sequence_desc>" + getSeqDesc() + "</sequence_desc>");
  str_vec.push_back("<proteoform_mass>" + std::to_string(getProteoformMass()) + "</proteoform_mass>");
  str_vec.push_back("<score>" + std::to_string(getScore()) + "</score>");

  if (n_trunc_shifts_.size() == 0) {
    str_vec.push_back("<n_trunc_shift_list/>");
  } else {
    str_vec.push_back("<n_trunc_shift_list>");
    for (size_t i = 0; i < n_trunc_shifts_.size(); i++) {
      str_vec.push_back("<shift>" + std::to_string(n_trunc_shifts_[i]) + "</shift>");
    }
    str_vec.push_back("</n_trunc_shift_list>");
  }

  if (c_trunc_shifts_.size() == 0) {
    str_vec.push_back("<c_trunc_shift_list/>");
  } else {
    str_vec.push_back("<c_trunc_shift_list>");
    for (size_t j = 0; j < c_trunc_shifts_.size(); j++) {
      str_vec.push_back("<shift>" + std::to_string(c_trunc_shifts_[j]) + "</shift>");
    }
    str_vec.push_back("</c_trunc_shift_list>");
  }

  str_vec.push_back("</simple_prsm>");
  return str_vec;
}

void SimplePrsm::setCTruncShifts(const std::vector<double> &c_term_shifts) {
  for (size_t i = 0; i < c_term_shifts.size(); i++) {
    double shift = prec_mass_ - (prot_mass_ + c_term_shifts[i]);
    c_trunc_shifts_.push_back(shift);
  }
}

} /* namespace prot */
