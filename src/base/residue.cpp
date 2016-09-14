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


#include <base/logger.hpp>

#include "base/acid_base.hpp"
#include "base/ptm_base.hpp"
#include "base/residue.hpp"
#include "base/file_util.hpp"
#include "base/string_util.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

Residue::Residue(AcidPtr acid_ptr, PtmPtr ptm_ptr):
    acid_ptr_(acid_ptr),
    ptm_ptr_(ptm_ptr) {
      mass_ = acid_ptr_->getMonoMass() + ptm_ptr_->getMonoMass();
    }

Residue::Residue(const std::string &acid_name, 
                 const std::string &ptm_abbr_name) {
  acid_ptr_ = AcidBase::getAcidPtrByName(acid_name);
  ptm_ptr_ = PtmBase::getPtmPtrByAbbrName(ptm_abbr_name);
  mass_ = acid_ptr_->getMonoMass() + ptm_ptr_->getMonoMass();
}

Residue::Residue(xercesc::DOMElement* element) { 
  std::string acid_element_name = Acid::getXmlElementName();
  xercesc::DOMElement* acid_element 
      = XmlDomUtil::getChildElement(element, acid_element_name.c_str(), 0);
  acid_ptr_ = AcidBase::getAcidPtrFromXml(acid_element);
  std::string ptm_element_name = Ptm::getXmlElementName();
  xercesc::DOMElement* ptm_element 
      = XmlDomUtil::getChildElement(element, ptm_element_name.c_str(), 0);
  ptm_ptr_ = PtmBase::getPtmPtrFromXml(ptm_element);  
  mass_ = acid_ptr_->getMonoMass() + ptm_ptr_->getMonoMass();
}

std::string Residue::toString(const std::string &delim_bgn, 
                              const std::string &delim_end) {
  if (PtmBase::isEmptyPtmPtr(ptm_ptr_)) {
    return acid_ptr_->getOneLetter();
  } else {
    return acid_ptr_->getOneLetter() + delim_bgn + ptm_ptr_->getAbbrName()
        + delim_end;
  }
}

void Residue::appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent,
                        const std::string &element_name){
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  std::string str = StringUtil::convertToString(mass_);
  acid_ptr_->appendNameToXml(xml_doc,element);
  ptm_ptr_->appendAbbrNameToXml(xml_doc,element);
  parent->appendChild(element);
}

void Residue::appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
  std::string element_name = Residue::getXmlElementName();
  appendXml(xml_doc, parent, element_name);
}

}
