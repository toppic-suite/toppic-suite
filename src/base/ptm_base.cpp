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

#include "base/ptm_base.hpp"
#include "base/logger.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

PtmPtrVec PtmBase::ptm_ptr_vec_;

PtmPtr PtmBase::empty_ptm_ptr_;

PtmPtr PtmBase::acetylation_ptr_;
PtmPtr PtmBase::c57_ptr_;
PtmPtr PtmBase::c58_ptr_;

void PtmBase::initBase(const std::string &file_name) {
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    std::string element_name = Ptm::getXmlElementName();
    int ptm_num = XmlDomUtil::getChildCount(parent, element_name.c_str());
    for (int i = 0; i < ptm_num; i++) {
      xercesc::DOMElement* element = XmlDomUtil::getChildElement(parent, element_name.c_str(), i);
      PtmPtr ptm_ptr(new Ptm(element));
      ptm_ptr_vec_.push_back(ptm_ptr);
      // check empty ptr
      if (ptm_ptr->getMonoMass() == 0.0) {
        empty_ptm_ptr_ = ptm_ptr;
      }
      // check acetylation
      if (ptm_ptr->getAbbrName() == PtmBase::getAcetylationAbbrName()) {
        acetylation_ptr_ = ptm_ptr;
      }
      if (ptm_ptr->getAbbrName() == PtmBase::getC57AbbrName()) {
        c57_ptr_ = ptm_ptr;
      }
      if (ptm_ptr->getAbbrName() == PtmBase::getC58AbbrName()) {
        c58_ptr_ = ptm_ptr;
      }
    }
  }
  if (empty_ptm_ptr_ == nullptr || acetylation_ptr_ == nullptr 
      || c57_ptr_ == nullptr || c58_ptr_ == nullptr) {
    LOG_WARN("ptm missing!");
  }
  std::sort(ptm_ptr_vec_.begin(), ptm_ptr_vec_.end(), Ptm::cmpMassInc);
}

// Returns a PTM based on the abbreviation name. Returns null if the
// abbreviation name does not exist.
PtmPtr PtmBase::getPtmPtrByAbbrName(const std::string &abbr_name) {
  for (size_t i = 0; i < ptm_ptr_vec_.size(); i++) {
    std::string n = ptm_ptr_vec_[i]->getAbbrName();
    if (n == abbr_name) {
      return ptm_ptr_vec_[i];
    }
  }
  return PtmPtr(nullptr);
}

PtmPtr PtmBase::getPtmPtr(PtmPtr p) {
  for (size_t i = 0; i < ptm_ptr_vec_.size(); i++) {
    if (ptm_ptr_vec_[i]->isSame(p)) {
      return ptm_ptr_vec_[i];
    }
  }
  ptm_ptr_vec_.push_back(p);
  return p;
}

// Checks if the list contains an amino acid with the specific name.
bool PtmBase::containsAbbrName(const std::string &abbr_name) {
  return getPtmPtrByAbbrName(abbr_name).get() != nullptr;
}

PtmPtr PtmBase::getPtmPtrFromXml(xercesc::DOMElement * element) {
  std::string abbr_name = Ptm::getAbbrNameFromXml(element);
  PtmPtr ptm_ptr = getPtmPtrByAbbrName(abbr_name);
  return ptm_ptr;
}

}

