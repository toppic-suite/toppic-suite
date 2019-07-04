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

#include "xercesc/framework/MemBufInputSource.hpp"

#include "common/util/logger.hpp"
#include "common/xml/xml_dom_document.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "common/base/amino_acid_base.hpp"
#include "common/base/amino_acid_data.hpp"

namespace toppic {

AminoAcidPtrVec AminoAcidBase::amino_acid_ptr_vec_;

AminoAcidPtr AminoAcidBase::empty_amino_acid_ptr_;

std::unordered_map<std::string, AminoAcidPtr> AminoAcidBase::amino_acid_one_letter_map_;

std::unordered_map<std::string, AminoAcidPtr> AminoAcidBase::amino_acid_three_letter_map_;

std::unordered_map<std::string, AminoAcidPtr> AminoAcidBase::amino_acid_name_map_;

void AminoAcidBase::initBase() {
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (!parser) {
    LOG_ERROR("Error in parsing amino acid data!");
    exit(EXIT_FAILURE);
  }

  xercesc::MemBufInputSource mem_str((const XMLByte*)amino_acid_base_data.c_str(), 
                                     amino_acid_base_data.length(), 
                                     "amino_acid_data");
  XmlDOMDocument doc(parser,mem_str);
  XmlDOMElement* parent = doc.getDocumentElement();
  std::string element_name = AminoAcid::getXmlElementName();
  int acid_num = xml_dom_util::getChildCount(parent, element_name.c_str());
  LOG_DEBUG("acid num " << acid_num);
  for (int i = 0; i < acid_num; i++) {
    XmlDOMElement* element = xml_dom_util::getChildElement(parent, element_name.c_str(), i);
    AminoAcidPtr ptr = std::make_shared<AminoAcid>(element);
    amino_acid_ptr_vec_.push_back(ptr);

    amino_acid_one_letter_map_[ptr->getOneLetter()]     = ptr;
    amino_acid_three_letter_map_[ptr->getThreeLetter()] = ptr;
    amino_acid_name_map_[ptr->getName()]                = ptr;

    // check if it is an empty acid
    if (ptr->getMonoMass() == 0.0) {
      empty_amino_acid_ptr_ = ptr;
    }
  }
}


AminoAcidPtr AminoAcidBase::getAminoAcidPtrByName(const std::string &name) {
  return amino_acid_name_map_[name];
}

AminoAcidPtr AminoAcidBase::getAminoAcidPtrByOneLetter(const std::string &one_letter) {
  return amino_acid_one_letter_map_[one_letter];
}

AminoAcidPtr AminoAcidBase::getAminoAcidPtrByThreeLetter(const std::string &three_letter) {
  return amino_acid_three_letter_map_[three_letter];
}

bool AminoAcidBase::containsName(const std::string &name) {
  return getAminoAcidPtrByName(name).get() != nullptr;
}

bool AminoAcidBase::containsOneLetter(const std::string &one_letter) {
  return getAminoAcidPtrByOneLetter(one_letter).get() != nullptr;
}

bool AminoAcidBase::containsThreeLetter(const std::string &three_letter) {
  return getAminoAcidPtrByThreeLetter(three_letter).get() != nullptr;
}

AminoAcidPtr AminoAcidBase::getAminoAcidPtrFromXml(XmlDOMElement * element) {
  std::string name = AminoAcid::getNameFromXml(element);
  AminoAcidPtr acid_ptr = getAminoAcidPtrByName(name);
  return acid_ptr;
}

}  // namespace toppic
