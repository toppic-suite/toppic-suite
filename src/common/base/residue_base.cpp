//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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

#include <string>

#include "common/util/logger.hpp"
#include "common/xml/xml_dom_document.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "common/base/amino_acid_base.hpp"
#include "common/base/ptm_base.hpp"
#include "common/base/residue_data.hpp"
#include "common/base/residue_base.hpp"

namespace toppic {

ResiduePtrVec ResidueBase::residue_ptr_vec_;
ResiduePtr ResidueBase::empty_residue_ptr_;

void ResidueBase::initBase() {
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (!parser) {
    LOG_ERROR("Error in parsing residue data!");
    exit(EXIT_FAILURE);
  }

  xercesc::MemBufInputSource mem_str((const XMLByte*)residue_base_data.c_str(), 
                                     residue_base_data.length(), 
                                     "residue_data");
  XmlDOMDocument doc(parser, mem_str);
  XmlDOMElement* parent = doc.getDocumentElement();
  std::string element_name = Residue::getXmlElementName();
  int residue_num = xml_dom_util::getChildCount(parent, element_name.c_str());
  LOG_DEBUG("residue num " << residue_num);
  for (int i = 0; i < residue_num; i++) {
    XmlDOMElement* element
      = xml_dom_util::getChildElement(parent, element_name.c_str(), i);
    ResiduePtr residue_ptr = std::make_shared<Residue>(element);
    if (residue_ptr->getAminoAcidPtr() == AminoAcidBase::getEmptyAminoAcidPtr()
        && residue_ptr->getPtmPtr() == PtmBase::getEmptyPtmPtr()) {
      empty_residue_ptr_ = residue_ptr;
    }
    residue_ptr_vec_.push_back(residue_ptr);
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

ResiduePtr ResidueBase::getBaseResiduePtr(AminoAcidPtr acid_ptr, PtmPtr ptm_ptr) {
  ResiduePtr residue_ptr = std::make_shared<Residue>(acid_ptr, ptm_ptr);
  return getBaseResiduePtr(residue_ptr);
}

ResiduePtr ResidueBase::getBaseResiduePtr(AminoAcidPtr acid_ptr) {
  ResiduePtr residue_ptr = std::make_shared<Residue>(acid_ptr, PtmBase::getEmptyPtmPtr());
  return getBaseResiduePtr(residue_ptr);
}

ResiduePtr ResidueBase::getResiduePtrFromXml(XmlDOMElement * element) {
  ResiduePtr ptr = std::make_shared<Residue>(element);
  return getBaseResiduePtr(ptr);
}

}  // namespace toppic
