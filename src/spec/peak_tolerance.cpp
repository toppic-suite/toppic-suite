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


#include "base/string_util.hpp"
#include "base/xml_dom_util.hpp"
#include "spec/peak_tolerance.hpp"


namespace prot {

PeakTolerance::PeakTolerance(double ppo, bool use_min_tolerance,
                             double min_tolerance):
    ppo_(ppo), 
    use_min_tolerance_(use_min_tolerance),
    min_tolerance_(min_tolerance) {
    }

PeakTolerance::PeakTolerance(xercesc::DOMElement* element){
  ppo_ = XmlDomUtil::getDoubleChildValue(element,"ppo",0);
  use_min_tolerance_ = XmlDomUtil::getDoubleChildValue(element,"use_min_tolerance",0);
  min_tolerance_ = XmlDomUtil::getDoubleChildValue(element,"min_tolerance",0);
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
  std::string str = StringUtil::convertToString(ppo_);
  xml_doc->addElement(element, "ppo", str.c_str());
  str = StringUtil::convertToString(use_min_tolerance_);
  xml_doc->addElement(element, "use_min_tolerance", str.c_str());
  str = StringUtil::convertToString(min_tolerance_);
  xml_doc->addElement(element, "min_tolerance", str.c_str());
  parent->appendChild(element);
}

} /* namespace prot */
