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


#include <cmath>
#include <limits>
#include <string>

#include "base/extreme_value.hpp"
#include "base/string_util.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

ExtremeValue::ExtremeValue(double one_prot_prob, double test_num,
                           double adjust_factor) {
  one_prot_prob_ = one_prot_prob;
  test_num_ = test_num;
  adjust_factor_ = adjust_factor;
  init();
}

void ExtremeValue::setOneProtProb(double one_prot_prob) {
  one_prot_prob_ = one_prot_prob;
  init();
}

ExtremeValue::ExtremeValue(xercesc::DOMElement* element) {
  one_prot_prob_ = XmlDomUtil::getDoubleChildValue(element, "one_protein_probability", 0);
  test_num_ = XmlDomUtil::getDoubleChildValue(element, "test_number", 0);
  adjust_factor_ = XmlDomUtil::getDoubleChildValue(element, "adjust_factor", 0);
  init();
}

void ExtremeValue::init() {
  e_value_ = one_prot_prob_ * test_num_ * adjust_factor_;
  if (one_prot_prob_ > 1 || test_num_ == std::numeric_limits<double>::max()) {
    p_value_  = 1.0;
  } else {
    double n = test_num_ * adjust_factor_;
    // approximation of 1 - (1- one_prot_prob)^n
    p_value_ =  n * one_prot_prob_
        - (n * (n - 1)) / 2 * one_prot_prob_ * one_prot_prob_
        + (n * (n - 1) * (n - 2)) / 6 * std::pow(one_prot_prob_, 3);
    if (p_value_ > 1.0) {
      p_value_ = 1.0;
    }
  }
}

void ExtremeValue::appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent) {
  std::string element_name = getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  std::string str = StringUtil::convertToString(one_prot_prob_);
  xml_doc->addElement(element, "one_protein_probability", str.c_str());
  str = StringUtil::convertToString(test_num_);
  xml_doc->addElement(element, "test_number", str.c_str());
  str = StringUtil::convertToString(adjust_factor_);
  xml_doc->addElement(element, "adjust_factor", str.c_str());
  str = StringUtil::convertToString(p_value_);
  xml_doc->addElement(element, "p_value", str.c_str());
  str = StringUtil::convertToString(e_value_);
  xml_doc->addElement(element, "e_value", str.c_str());
  parent->appendChild(element);
}

ExtremeValuePtr ExtremeValue::getMaxEvaluePtr() {
  ExtremeValuePtr evalue_ptr
      = std::make_shared<ExtremeValue>(1.0, std::numeric_limits<double>::max(), 1.0);
  return evalue_ptr;
}

}  // namespace prot
