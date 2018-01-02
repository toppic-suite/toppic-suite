//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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
    int mod_num = xml_dom_util::getChildCount(parent, element_name.c_str());
    for (int i = 0; i < mod_num; i++) {
      xercesc::DOMElement* element = xml_dom_util::getChildElement(parent, element_name.c_str(), i);
      ModPtr mod_ptr = std::make_shared<Mod>(element);
      mod_ptr_vec_.push_back(mod_ptr);
      // check empty ptr
      if (mod_ptr->getOriResiduePtr() == ResidueBase::getEmptyResiduePtr()
          && mod_ptr->getModResiduePtr() ==ResidueBase::getEmptyResiduePtr()) {
        none_mod_ptr_ = mod_ptr;
      }
      if (mod_ptr->getModResiduePtr()->getAminoAcidPtr()->getOneLetter() == "C"
          && mod_ptr->getModResiduePtr()->getPtmPtr() == PtmBase::getPtmPtr_C57()) {
        c57_mod_ptr_ = mod_ptr;
      }
      if (mod_ptr->getModResiduePtr()->getAminoAcidPtr()->getOneLetter() == "C"
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
  ModPtr mod_ptr = std::make_shared<Mod>(ori_residue, mod_residue);
  return getBaseModPtr(mod_ptr);
}

ModPtr ModBase::getModPtrFromXml(xercesc::DOMElement * element) {
  ModPtr ptr = std::make_shared<Mod>(element);
  return getBaseModPtr(ptr);
}

}  // namespace prot

