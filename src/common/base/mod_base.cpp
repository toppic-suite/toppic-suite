//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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
#include "common/base/ptm_base.hpp"
#include "common/base/residue_base.hpp"
#include "common/base/mod_data.hpp"
#include "common/base/mod_base.hpp"

namespace toppic {

ModPtrVec ModBase::mod_ptr_vec_;
ModPtr ModBase::none_mod_ptr_;
ModPtr ModBase::c57_mod_ptr_;
ModPtr ModBase::c58_mod_ptr_;

void ModBase::initBase() {
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (!parser) {
    LOG_ERROR("Error in parsing modification data!");
    exit(EXIT_FAILURE);
  }

  xercesc::MemBufInputSource mem_str((const XMLByte*)mod_base_data.c_str(), 
                                     mod_base_data.length(), 
                                     "modification_data");
  XmlDOMDocument doc(parser, mem_str);
  XmlDOMElement* parent = doc.getDocumentElement();
  std::string element_name = Mod::getXmlElementName();
  int mod_num = xml_dom_util::getChildCount(parent, element_name.c_str());
  for (int i = 0; i < mod_num; i++) {
    XmlDOMElement* element = xml_dom_util::getChildElement(parent, element_name.c_str(), i);
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
    LOG_ERROR("Modification configuration file is incomplete!");
    exit(EXIT_FAILURE);
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

ModPtr ModBase::getModPtrFromXml(XmlDOMElement * element) {
  ModPtr ptr = std::make_shared<Mod>(element);
  return getBaseModPtr(ptr);
}

}  // namespace toppic

