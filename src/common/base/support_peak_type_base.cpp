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

#include <stdexcept>
#include <string>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "common/xml/xml_dom_parser.hpp"
#include "common/xml/xml_dom_document.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "common/base/support_peak_type_base.hpp"

namespace toppic {

SPTypePtrVec SPTypeBase::sp_type_ptr_vec_;
std::unordered_map<std::string, SPTypePtr> SPTypeBase::sp_type_name_map_;
std::unordered_map<int, SPTypePtr> SPTypeBase::sp_type_id_map_;

SPTypePtr SPTypeBase::sp_type_ptr_N_TERM_;

void SPTypeBase::initBase(const std::string &base_dir) {
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (!parser) {
    LOG_ERROR("Error in parsing support peak type data!");
    throw std::runtime_error("Error in parsing support peak type data!");
  }

  std::string sp_type_base_file_name = base_dir 
      + file_util::getFileSeparator() + "support_peak_type_base.xml";
  std::string sp_type_base_data = file_util::readFile(sp_type_base_file_name);
  xercesc::MemBufInputSource mem_str((const XMLByte*)sp_type_base_data.c_str(), 
                                     sp_type_base_data.length(), 
                                     "support_peak_type_data");
  XmlDOMDocument doc(parser, mem_str);
  XmlDOMElement* parent = doc.getDocumentElement();
  std::string element_name = SupportPeakType::getXmlElementName();
  int prm_peak_type_num = xml_dom_util::getChildCount(parent, element_name.c_str());
  for (int i = 0; i < prm_peak_type_num; i++) {
    XmlDOMElement* element = xml_dom_util::getChildElement(parent, element_name.c_str(), i);
    SPTypePtr sp_type_ptr = std::make_shared<SupportPeakType>(element);
    sp_type_ptr_vec_.push_back(sp_type_ptr);
    sp_type_name_map_[sp_type_ptr->getName()] = sp_type_ptr;
    sp_type_id_map_[sp_type_ptr->getId()] = sp_type_ptr;
    if (sp_type_ptr->getName() == getName_N_TERM()) {
      sp_type_ptr_N_TERM_ = sp_type_ptr;
    }
  }
}

SPTypePtr SPTypeBase::getSPTypePtrByName(const std::string &name) {
  auto it = sp_type_name_map_.find(name);
  if (it == sp_type_name_map_.end()) {
    LOG_WARN("Support peak type " << name << " cannot be found!");
    return nullptr;
  }
  return it->second;
}

SPTypePtr SPTypeBase::getSPTypePtrById(int id) {
  auto it = sp_type_id_map_.find(id);
  if (it == sp_type_id_map_.end()) {
    LOG_WARN("Support peak id " << id << " cannot be found!");
    return nullptr;
  }
  return it->second;
}

}  // namespace toppic
