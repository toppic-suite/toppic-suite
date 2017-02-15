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
#include "base/ptm_base.hpp"
#include "base/trunc_base.hpp"
#include "base/mod_base.hpp"
#include "base/prot_mod.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

ProtMod::ProtMod(const std::string &name, const std::string &type,
                 TruncPtr trunc_ptr, ModPtr mod_ptr): 
    name_(name),
    type_(type),
    trunc_ptr_(trunc_ptr),
    mod_ptr_(mod_ptr) {
      mod_pos_ = trunc_ptr->getTruncLen();
      prot_shift_ = trunc_ptr_->getShift() + mod_ptr_->getShift();
      pep_shift_ = mod_ptr_->getShift();
    }

ProtMod::ProtMod(xercesc::DOMElement* element) { 
  name_ = XmlDomUtil::getChildValue(element, "name", 0);
  type_ = XmlDomUtil::getChildValue(element, "type", 0);
  std::string trunc_element_name = Trunc::getXmlElementName();
  xercesc::DOMElement* trunc_element 
      = XmlDomUtil::getChildElement(element, trunc_element_name.c_str(), 0);
  trunc_ptr_ = TruncBase::getTruncPtrFromXml(trunc_element);
  std::string mod_element_name = Mod::getXmlElementName();
  xercesc::DOMElement* mod_element 
      = XmlDomUtil::getChildElement(element, mod_element_name.c_str(), 0);
  mod_ptr_= ModBase::getModPtrFromXml(mod_element); 
  mod_pos_ = trunc_ptr_->getTruncLen();
  prot_shift_ = trunc_ptr_->getShift() + mod_ptr_->getShift();
  pep_shift_ = mod_ptr_->getShift();
}

void ProtMod::appendNameToXml(XmlDOMDocument* xml_doc,
                              xercesc::DOMElement* parent){
  std::string element_name = ProtMod::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  xml_doc->addElement(element, "name", name_.c_str());
  parent->appendChild(element);
}

std::string ProtMod::getNameFromXml(xercesc::DOMElement * element) {
  std::string name = XmlDomUtil::getChildValue(element, "name", 0);
  return name;
}

bool ProtMod::isAcetylation() {
  if (mod_ptr_->getModResiduePtr()->getPtmPtr() == PtmBase::getPtmPtr_Acetylation()) {
    return true;
  }
  else {
    return false;
  }
}

}
