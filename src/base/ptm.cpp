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


#include <iostream>
#include <vector>
#include <memory>
#include <algorithm>

#include "base/ptm.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"
#include "base/xml_dom_util.hpp"
#include "base/string_util.hpp"

namespace prot {

Ptm::Ptm(const std::string &name, const std::string &abbr_name,
         double mono_mass, int unimod_id, 
         const std::string &n_term_residue_str,
         const std::string &c_term_residue_str, 
         const std::string &anywhere_residue_str): 
    name_(name),
    abbr_name_(abbr_name),
    mono_mass_(mono_mass),
    unimod_id_(unimod_id) {
    }

Ptm::Ptm(xercesc::DOMElement* element) { 
  name_ = XmlDomUtil::getChildValue(element, "name", 0);
  abbr_name_ = XmlDomUtil::getChildValue(element, "abbreviation", 0);
  mono_mass_ = XmlDomUtil::getDoubleChildValue(element, "mono_mass", 0);
  unimod_id_ = XmlDomUtil::getIntChildValue(element, "unimod", 0);
}

void Ptm::appendAbbrNameToXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent) {
  std::string element_name = Ptm::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  xml_doc->addElement(element, "abbreviation", abbr_name_.c_str());
  xml_doc->addElement(element, "unimod", StringUtil::convertToString(unimod_id_).c_str());
  parent->appendChild(element);
}

std::string Ptm::getAbbrNameFromXml(xercesc::DOMElement * element) {
  std::string abbr_name = XmlDomUtil::getChildValue(element, "abbreviation", 0);
  return abbr_name;
}

}

