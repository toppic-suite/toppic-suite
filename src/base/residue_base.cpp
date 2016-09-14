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


#include "base/logger.hpp"
#include "base/acid_base.hpp"
#include "base/ptm_base.hpp"
#include "base/residue_base.hpp"
#include "base/file_util.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

ResiduePtrVec ResidueBase::residue_ptr_vec_;
ResiduePtr ResidueBase::empty_residue_ptr_;

void ResidueBase::initBase(const std::string &file_name) {
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    std::string element_name = Residue::getXmlElementName();
    int residue_num = XmlDomUtil::getChildCount(parent, element_name.c_str());
    LOG_DEBUG( "residue num " << residue_num);
    for (int i = 0; i < residue_num; i++) {
      xercesc::DOMElement* element
          = XmlDomUtil::getChildElement(parent, element_name.c_str(), i);
      ResiduePtr residue_ptr(new Residue(element));
      if (residue_ptr->getAcidPtr() == AcidBase::getEmptyAcidPtr() 
          && residue_ptr->getPtmPtr() == PtmBase::getEmptyPtmPtr()) {
        empty_residue_ptr_ = residue_ptr;
      }
      residue_ptr_vec_.push_back(residue_ptr);
    }
  }
}

// use residue in residue_base to remove duplications and reduce memory usage
// add the residue to base if it is a new one
ResiduePtr ResidueBase::getBaseResiduePtr(ResiduePtr residue_ptr) {
  for (size_t i = 0; i < residue_ptr_vec_.size(); i++) {
    if (residue_ptr_vec_[i]->isSame(residue_ptr)) {
      return residue_ptr_vec_[i];
    }
  }
  residue_ptr_vec_.push_back(residue_ptr);
  return residue_ptr;
}

ResiduePtrVec ResidueBase::getBaseNonePtmResiduePtrVec() {
  ResiduePtrVec result;
  for (size_t i = 0; i < residue_ptr_vec_.size(); i++) {
    if (residue_ptr_vec_[i]->getPtmPtr() == PtmBase::getEmptyPtmPtr()) {
      result.push_back(residue_ptr_vec_[i]);
    }
  }
  return result;
}

ResiduePtr ResidueBase::getBaseResiduePtr(AcidPtr acid_ptr, PtmPtr ptm_ptr) {
  ResiduePtr residue_ptr = ResiduePtr(new Residue(acid_ptr, ptm_ptr));
  return getBaseResiduePtr(residue_ptr);
}

ResiduePtr ResidueBase::getBaseResiduePtr(AcidPtr acid_ptr) {
  ResiduePtr residue_ptr = ResiduePtr(new Residue(acid_ptr, PtmBase::getEmptyPtmPtr()));
  return getBaseResiduePtr(residue_ptr);
}

ResiduePtr ResidueBase::getResiduePtrFromXml(xercesc::DOMElement * element) {
  ResiduePtr ptr(new Residue(element));
  return getBaseResiduePtr(ptr);
}

}
