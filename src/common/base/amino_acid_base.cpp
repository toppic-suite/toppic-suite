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
#include "common/xml/xml_dom_document.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "common/base/amino_acid_base.hpp"

namespace toppic {

AminoAcidPtrVec AminoAcidBase::amino_acid_ptr_vec_;

AminoAcidPtr AminoAcidBase::empty_amino_acid_ptr_;

std::unordered_map<std::string, AminoAcidPtr> AminoAcidBase::amino_acid_one_letter_map_;

std::unordered_map<std::string, AminoAcidPtr> AminoAcidBase::amino_acid_three_letter_map_;

std::unordered_map<std::string, AminoAcidPtr> AminoAcidBase::amino_acid_name_map_;

void AminoAcidBase::initBase(const std::string &base_dir) {
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (!parser) {
    LOG_ERROR("Error in parsing amino acid data!");
    throw std::runtime_error("Error in parsing amino acid data!");
  }
  std::string amino_acid_base_file_name = base_dir 
      + file_util::getFileSeparator() + "amino_acid_base.xml";
  std::string amino_acid_base_data = file_util::readFile(amino_acid_base_file_name);

  xercesc::MemBufInputSource mem_str((const XMLByte*)amino_acid_base_data.c_str(), 
                                     amino_acid_base_data.length(), 
                                     "amino_acid_data");
  XmlDOMDocument doc(parser, mem_str);
  XmlDOMElement* parent = doc.getDocumentElement();
  std::string element_name = AminoAcid::getXmlElementName();
  int acid_num = xml_dom_util::getChildCount(parent, element_name.c_str());
  LOG_DEBUG("acid num " << acid_num);
  for (int i = 0; i < acid_num; i++) {
    XmlDOMElement* element 
        = xml_dom_util::getChildElement(parent, element_name.c_str(), i);
    AminoAcidPtr ptr = std::make_shared<AminoAcid>(element);
    amino_acid_ptr_vec_.push_back(ptr);

    amino_acid_one_letter_map_[ptr->getOneLetter()]     = ptr;
    amino_acid_three_letter_map_[ptr->getThreeLetter()] = ptr;
    amino_acid_name_map_[ptr->getName()]                = ptr;

    // check if it is an empty acid
    if (ptr->getName() == "None") {
      empty_amino_acid_ptr_ = ptr;
    }
  }
  if (empty_amino_acid_ptr_ == nullptr) {
    LOG_ERROR("Empty amino acid cannot be found in amino acid base data!");
    throw std::runtime_error("Empty amino acid cannot be found in amino acid base data!");
  }
}

AminoAcidPtr AminoAcidBase::getAminoAcidPtrByName(const std::string &name) {
  auto it = amino_acid_name_map_.find(name);
  if (it == amino_acid_name_map_.end()) {
    LOG_WARN("Amino acid " << name << " cannot be found!");
    return nullptr;
  }
  return it->second;
}

AminoAcidPtr AminoAcidBase::getAminoAcidPtrByOneLetter(const std::string &one_letter) {
  auto it = amino_acid_one_letter_map_.find(one_letter);
  if (it == amino_acid_one_letter_map_.end()) {
    LOG_WARN("Amino acid " << one_letter << " cannot be found!");
    return nullptr;
  }
  return it->second;
}

AminoAcidPtr
AminoAcidBase::getAminoAcidPtrByThreeLetter(const std::string &three_letter) {
  auto it = amino_acid_three_letter_map_.find(three_letter);
  if (it == amino_acid_three_letter_map_.end()) {
    LOG_WARN("Amino acid " << three_letter << " cannot be found!");
    return nullptr;
  }
  return it->second;
}

bool AminoAcidBase::containsName(const std::string &name) {
  return getAminoAcidPtrByName(name) != nullptr;
}

bool AminoAcidBase::containsOneLetter(const std::string &one_letter) {
  return getAminoAcidPtrByOneLetter(one_letter) != nullptr;
}

bool AminoAcidBase::containsThreeLetter(const std::string &three_letter) {
  return getAminoAcidPtrByThreeLetter(three_letter) != nullptr;
}

AminoAcidPtr AminoAcidBase::getAminoAcidPtrFromXml(XmlDOMElement * element) {
  std::string name = AminoAcid::getNameFromXml(element);
  return getAminoAcidPtrByName(name);
}

}  // namespace toppic
