//Copyright (c) 2014 - 2025, The Trustees of Indiana University, Tulane University.
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
#include "common/util/file_util.hpp"
#include "common/xml/xml_dom_document.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "common/base/ptm_base.hpp"
#include "common/base/residue_base.hpp"
#include "common/base/n_term_mod_base.hpp"

namespace toppic {

ModPtrVec NTermModBase::mod_ptr_vec_;
ModPtr NTermModBase::none_mod_ptr_;

void NTermModBase::initBase(const std::string &base_dir) {
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (!parser) {
    LOG_ERROR("Error in parsing modification data!");
    exit(EXIT_FAILURE);
  }

  std::string mod_base_file_name = base_dir 
      + file_util::getFileSeparator() + "n_term_mod_base.xml";
  std::string mod_base_data = file_util::readFile(mod_base_file_name);
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
    if (mod_ptr->getOriResiduePtr() == ResidueBase::getEmptyResiduePtr() &&
        mod_ptr->getModResiduePtr() == ResidueBase::getEmptyResiduePtr()) {
      none_mod_ptr_ = mod_ptr;
    }
  }
  if (none_mod_ptr_ == nullptr) {
    LOG_ERROR("N-term modification configuration file is incomplete!");
    exit(EXIT_FAILURE);
  }
}

ModPtr NTermModBase::getBaseNTermModPtr(ModPtr n_term_mod_ptr) {
  for (size_t i = 0; i < mod_ptr_vec_.size(); i++) {
    if (mod_ptr_vec_[i]->isSame(n_term_mod_ptr)) {
      return mod_ptr_vec_[i];
    }
  }
  mod_ptr_vec_.push_back(n_term_mod_ptr);
  return n_term_mod_ptr;
}


ModPtr NTermModBase::getNTermModPtrFromXml(XmlDOMElement * element) {
  ModPtr ptr = std::make_shared<Mod>(element);
  return getBaseNTermModPtr(ptr);
}

}  // namespace toppic

