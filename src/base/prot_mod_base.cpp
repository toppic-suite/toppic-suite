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
#include "base/prot_mod_base.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

ProtModPtrVec ProtModBase::prot_mod_ptr_vec_;
ProtModPtr ProtModBase::prot_mod_ptr_NONE_;
ProtModPtr ProtModBase::prot_mod_ptr_M_ACETYLATION_;
//ProtModPtr ProtModBase::prot_mod_ptr_NME_;
//ProtModPtr ProtModBase::prot_mod_ptr_NME_ACETYLATION_;

void ProtModBase::initBase(const std::string &file_name) {
  prot::XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    prot::XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    std::string element_name = ProtMod::getXmlElementName();
    int mod_num = XmlDomUtil::getChildCount(parent, element_name.c_str());
    for (int i = 0; i < mod_num; i++) {
      xercesc::DOMElement* element = XmlDomUtil::getChildElement(parent, element_name.c_str(), i);
      ProtModPtr prot_mod_ptr(new ProtMod(element));
      //LOG_DEBUG("ptm index " << i << " shift  " << prot_mod_ptr->getProtShift());
      prot_mod_ptr_vec_.push_back(prot_mod_ptr);
      if (prot_mod_ptr->getName() == getName_NONE()) {
        prot_mod_ptr_NONE_ = prot_mod_ptr;
      }
      if (prot_mod_ptr->getName() == getName_M_ACETYLATION()) {
        prot_mod_ptr_M_ACETYLATION_ = prot_mod_ptr;
      }
      /*
      if (prot_mod_ptr->getName() == getName_NME()) {
        prot_mod_ptr_NME_ = prot_mod_ptr;
      }
      if (prot_mod_ptr->getName() == getName_NME_ACETYLATION()) {
        prot_mod_ptr_NME_ACETYLATION_ = prot_mod_ptr;
      }
      */
    }
  }
}

ProtModPtr ProtModBase::getProtModPtrByName(const std::string &name) {
  for (size_t i = 0; i < prot_mod_ptr_vec_.size(); i++) {
    std::string n = prot_mod_ptr_vec_[i]->getName();
    if (n == name) {
      return prot_mod_ptr_vec_[i];
    }
  }
  LOG_WARN("prot mod nullptr");
  return ProtModPtr(nullptr);
}

ProtModPtrVec ProtModBase::getProtModPtrByType(const std::string &type) {
  ProtModPtrVec prot_mods;
  for (size_t i = 0; i < prot_mod_ptr_vec_.size(); i++) {
    std::string t = prot_mod_ptr_vec_[i]->getType();
    if (t == type) {
      prot_mods.push_back(prot_mod_ptr_vec_[i]);
    }
  }
  return prot_mods;
}

ProtModPtr ProtModBase::getProtModPtrFromXml(xercesc::DOMElement * element) {
  std::string name = ProtMod::getNameFromXml(element);
  ProtModPtr prot_mod_ptr = getProtModPtrByName(name);
  return prot_mod_ptr;
}

}
