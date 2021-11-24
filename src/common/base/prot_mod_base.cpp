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
#include "common/xml/xml_dom_parser.hpp"
#include "common/xml/xml_dom_document.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "common/base/prot_mod_data.hpp"
#include "common/base/prot_mod_base.hpp"

namespace toppic {

ProtModPtrVec ProtModBase::prot_mod_ptr_vec_;
ProtModPtr ProtModBase::prot_mod_ptr_NONE_;
ProtModPtr ProtModBase::prot_mod_ptr_M_ACETYLATION_;
ProtModPtr ProtModBase::prot_mod_ptr_NME_;

void ProtModBase::initBase() {
  toppic::XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (!parser) {
    LOG_ERROR("Error in parsing protein modification data!");
    exit(EXIT_FAILURE);
  }

  xercesc::MemBufInputSource mem_str((const XMLByte*)prot_mod_base_data.c_str(), 
                                     prot_mod_base_data.length(), 
                                     "prot_mod_data");
  XmlDOMDocument doc(parser, mem_str);
  XmlDOMElement* parent = doc.getDocumentElement();
  std::string element_name = ProtMod::getXmlElementName();
  int mod_num = xml_dom_util::getChildCount(parent, element_name.c_str());
  for (int i = 0; i < mod_num; i++) {
    XmlDOMElement* element = xml_dom_util::getChildElement(parent, element_name.c_str(), i);
    ProtModPtr prot_mod_ptr = std::make_shared<ProtMod>(element);
    prot_mod_ptr_vec_.push_back(prot_mod_ptr);
    if (prot_mod_ptr->getType() == getType_NONE()) {
      prot_mod_ptr_NONE_ = prot_mod_ptr;
    }
    if (prot_mod_ptr->getType() == getType_M_ACETYLATION()) {
      prot_mod_ptr_M_ACETYLATION_ = prot_mod_ptr;
    }
    if (prot_mod_ptr->getType() == getType_NME()) {
      prot_mod_ptr_NME_ = prot_mod_ptr;
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
  LOG_WARN("Protein modification " << name << " cannot be found!");
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

ProtModPtr ProtModBase::getProtModPtrFromXml(XmlDOMElement * element) {
  std::string name = ProtMod::getNameFromXml(element);
  ProtModPtr prot_mod_ptr = getProtModPtrByName(name);
  return prot_mod_ptr;
}

}  // namespace toppic
