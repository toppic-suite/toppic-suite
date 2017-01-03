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

#include <string>

#include "base/logger.hpp"
#include "base/acid_base.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

AcidPtrVec AcidBase::acid_ptr_vec_;
AcidPtr AcidBase::empty_acid_ptr_;

void AcidBase::initBase(const std::string &file_name) {
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    std::string element_name = Acid::getXmlElementName();
    int acid_num = XmlDomUtil::getChildCount(parent, element_name.c_str());
    LOG_DEBUG("acid num " << acid_num);
    for (int i = 0; i < acid_num; i++) {
      xercesc::DOMElement* element
          = XmlDomUtil::getChildElement(parent, element_name.c_str(), i);
      AcidPtr ptr(new Acid(element));
      acid_ptr_vec_.push_back(ptr);
      // check if it is an empty acid
      if (ptr->getMonoMass() == 0.0) {
        empty_acid_ptr_ = ptr;
      }
    }
  }
}

AcidPtr AcidBase::getAcidPtrByName(const std::string &name) {
  for (size_t i = 0; i < acid_ptr_vec_.size(); i++) {
    std::string n = acid_ptr_vec_[i]->getName();
    if (n == name) {
      return acid_ptr_vec_[i];
    }
  }
  LOG_DEBUG("Acid not found: " + name);
  return AcidPtr(nullptr);
}

AcidPtr AcidBase::getAcidPtrByOneLetter(const std::string &one_letter) {
  for (size_t i = 0; i < acid_ptr_vec_.size(); i++) {
    std::string l = acid_ptr_vec_[i]->getOneLetter();
    if (l == one_letter)  {
      return acid_ptr_vec_[i];
    }
  }
  LOG_DEBUG("Acid not found " + one_letter);
  return AcidPtr(nullptr);
}

AcidPtr AcidBase::getAcidPtrByThreeLetter(const std::string &three_letter) {
  for (size_t i = 0; i < acid_ptr_vec_.size(); i++) {
    std::string l = acid_ptr_vec_[i]->getThreeLetter();
    if (l == three_letter) {
      return acid_ptr_vec_[i];
    }
  }
  LOG_DEBUG("Acid not found " + three_letter);
  return AcidPtr(nullptr);
}

bool AcidBase::containsName(const std::string &name) {
  return getAcidPtrByName(name).get() != nullptr;
}

bool AcidBase::containsOneLetter(const std::string &one_letter) {
  return getAcidPtrByOneLetter(one_letter).get() != nullptr;
}

bool AcidBase::containsThreeLetter(const std::string &three_letter) {
  return getAcidPtrByThreeLetter(three_letter).get() != nullptr;
}

AcidPtr AcidBase::getAcidPtrFromXml(xercesc::DOMElement * element) {
  std::string name = Acid::getNameFromXml(element);
  AcidPtr acid_ptr = getAcidPtrByName(name);
  return acid_ptr;
}

}  // namespace prot
