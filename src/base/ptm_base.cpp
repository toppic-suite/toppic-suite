//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.


#include <iostream>
#include <vector>
#include <memory>
#include <algorithm>
#include <string>

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
      PtmPtr ptm_ptr = std::make_shared<Ptm>(element);
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

}  // namespace prot

