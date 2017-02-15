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


#include "base/string_util.hpp"
#include "base/xml_dom_util.hpp"
#include "spec/deconv_peak.hpp"

namespace prot {

DeconvPeak::DeconvPeak (int id, double mono_mass, 
                        double intensity, int charge):
    Peak (mono_mass, intensity),
    id_(id),
    charge_(charge),
    score_(1.0) {
    }

DeconvPeak::DeconvPeak(xercesc::DOMElement* element):
    Peak (XmlDomUtil::getDoubleChildValue(element,"position",0), 
          XmlDomUtil::getDoubleChildValue(element,"intensity",0)) {
      id_ = XmlDomUtil::getIntChildValue(element,"id",0);
      charge_ = XmlDomUtil::getIntChildValue(element,"charge",0);
      score_ = XmlDomUtil::getDoubleChildValue(element,"score",0);
    }

void DeconvPeak::appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
  std::string element_name = DeconvPeak::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  std::string str = StringUtil::convertToString(getPosition());
  xml_doc->addElement(element, "position", str.c_str());
  str = StringUtil::convertToString(getIntensity());
  xml_doc->addElement(element, "intensity", str.c_str());
  str = StringUtil::convertToString(id_);
  xml_doc->addElement(element, "id", str.c_str());
  str = StringUtil::convertToString(charge_);
  xml_doc->addElement(element, "charge", str.c_str());
  str = StringUtil::convertToString(score_);
  xml_doc->addElement(element, "score", str.c_str());
  parent->appendChild(element);
}

}

