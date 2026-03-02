//Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane University.
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

#include <stdexcept>
#include <string>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "common/xml/xml_dom_parser.hpp"
#include "common/xml/xml_dom_document.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "common/base/prot_mod_base.hpp"

namespace toppic {

ProtModPtrVec ProtModBase::prot_mod_ptr_vec_;
std::unordered_map<std::string, ProtModPtr> ProtModBase::prot_mod_name_map_;
ProtModPtr ProtModBase::prot_mod_ptr_NONE_;
ProtModPtr ProtModBase::prot_mod_ptr_M_ACETYLATION_;
ProtModPtr ProtModBase::prot_mod_ptr_NME_;

void ProtModBase::initBase(const std::string &base_dir) {
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (!parser) {
    LOG_ERROR("Error in parsing protein modification data!");
    throw std::runtime_error("Error in parsing protein modification data!");
  }

  std::string prot_mod_base_file_name = base_dir 
      + file_util::getFileSeparator() + "prot_mod_base.xml";
  std::string prot_mod_base_data = file_util::readFile(prot_mod_base_file_name);
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
    prot_mod_name_map_[prot_mod_ptr->getName()] = prot_mod_ptr;
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
  auto it = prot_mod_name_map_.find(name);
  if (it == prot_mod_name_map_.end()) {
    LOG_ERROR("Protein modification " << name << " cannot be found!");
    return nullptr;
  }
  return it->second;
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
  return getProtModPtrByName(name);
}

}  // namespace toppic
