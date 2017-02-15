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


#include "base/logger.hpp"
#include "base/mod.hpp"
#include "base/residue_base.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

Mod::Mod(ResiduePtr ori_residue_ptr, ResiduePtr mod_residue_ptr):
    ori_residue_ptr_(ori_residue_ptr),
    mod_residue_ptr_(mod_residue_ptr) {
    }

Mod::Mod(xercesc::DOMElement* element) { 
  xercesc::DOMElement* ori_residue_element 
      = XmlDomUtil::getChildElement(element, "ori_residue", 0);
  ori_residue_ptr_ = ResidueBase::getResiduePtrFromXml(ori_residue_element);
  xercesc::DOMElement* mod_residue_element 
      = XmlDomUtil::getChildElement(element, "mod_residue", 0);
  mod_residue_ptr_ = ResidueBase::getResiduePtrFromXml(mod_residue_element);
}

void Mod::appendToXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent){
  std::string element_name = Mod::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  ori_residue_ptr_->appendXml(xml_doc, element, "ori_residue");
  mod_residue_ptr_->appendXml(xml_doc, element, "mod_residue");
  parent->appendChild(element);
}

}
