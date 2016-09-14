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
#include <vector>
#include <memory>
#include <algorithm>

#include "base/ptm_base.hpp"
#include "base/mod_base.hpp"
#include "base/logger.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

ModPtrVec ModBase::mod_ptr_vec_;
ModPtr ModBase::none_mod_ptr_;
ModPtr ModBase::c57_mod_ptr_;
ModPtr ModBase::c58_mod_ptr_;

void ModBase::initBase(const std::string &file_name) {
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    std::string element_name = Mod::getXmlElementName();
    int mod_num = XmlDomUtil::getChildCount(parent, element_name.c_str());
    for (int i = 0; i < mod_num; i++) {
      xercesc::DOMElement* element = XmlDomUtil::getChildElement(parent, element_name.c_str(), i);
      ModPtr mod_ptr(new Mod(element));
      mod_ptr_vec_.push_back(mod_ptr);
      // check empty ptr
      if (mod_ptr->getOriResiduePtr() == ResidueBase::getEmptyResiduePtr() 
          && mod_ptr->getModResiduePtr() ==ResidueBase::getEmptyResiduePtr()) {
        none_mod_ptr_ = mod_ptr;
      }
      if (mod_ptr->getModResiduePtr()->getAcidPtr()->getOneLetter() == "C"
          && mod_ptr->getModResiduePtr()->getPtmPtr() == PtmBase::getPtmPtr_C57()) {
        c57_mod_ptr_ = mod_ptr;
      }
      if (mod_ptr->getModResiduePtr()->getAcidPtr()->getOneLetter() == "C"
          && mod_ptr->getModResiduePtr()->getPtmPtr() == PtmBase::getPtmPtr_C58()) {
        c58_mod_ptr_ = mod_ptr;
      }
    }
    if (none_mod_ptr_ == nullptr || c57_mod_ptr_ == nullptr || c58_mod_ptr_ == nullptr) {
      LOG_WARN("mod missing!");
    }
  }
}

ModPtr ModBase::getBaseModPtr(ModPtr mod_ptr) {
  for (size_t i = 0; i < mod_ptr_vec_.size(); i++) {
    if (mod_ptr_vec_[i]->isSame(mod_ptr)) {
      return mod_ptr_vec_[i];
    }
  }
  mod_ptr_vec_.push_back(mod_ptr);
  return mod_ptr;
}

ModPtr ModBase::getBaseModPtr(ResiduePtr ori_residue, ResiduePtr mod_residue) {
  ModPtr mod_ptr(new Mod(ori_residue, mod_residue));
  return getBaseModPtr(mod_ptr);
}

ModPtr ModBase::getModPtrFromXml(xercesc::DOMElement * element) {
  ModPtr ptr(new Mod(element));
  return getBaseModPtr(ptr);
}
}

