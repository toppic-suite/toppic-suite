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

#include <string>
#include "base/logger.hpp"
#include "base/prot_mod_base.hpp"
#include "base/xml_dom_util.hpp"

namespace toppic {

ProtModPtrVec ProtModBase::prot_mod_ptr_vec_;
ProtModPtr ProtModBase::prot_mod_ptr_NONE_;
ProtModPtr ProtModBase::prot_mod_ptr_M_ACETYLATION_;
//  ProtModPtr ProtModBase::prot_mod_ptr_NME_;
//  ProtModPtr ProtModBase::prot_mod_ptr_NME_ACETYLATION_;

void ProtModBase::initBase(const std::string &file_name) {
  toppic::XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    toppic::XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    std::string element_name = ProtMod::getXmlElementName();
    int mod_num = xml_dom_util::getChildCount(parent, element_name.c_str());
    for (int i = 0; i < mod_num; i++) {
      xercesc::DOMElement* element = xml_dom_util::getChildElement(parent, element_name.c_str(), i);
      ProtModPtr prot_mod_ptr = std::make_shared<ProtMod>(element);
      //  LOG_DEBUG("ptm index " << i << " shift  " << prot_mod_ptr->getProtShift());
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

}  // namespace toppic
