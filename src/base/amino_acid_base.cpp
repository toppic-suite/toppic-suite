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
#include "base/amino_acid_base.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

AminoAcidPtrVec AminoAcidBase::amino_acid_ptr_vec_;
AminoAcidPtr AminoAcidBase::empty_amino_acid_ptr_;

void AminoAcidBase::initBase(const std::string &file_name) {
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    std::string element_name = AminoAcid::getXmlElementName();
    int acid_num = xml_dom_util::getChildCount(parent, element_name.c_str());
    LOG_DEBUG("acid num " << acid_num);
    for (int i = 0; i < acid_num; i++) {
      xercesc::DOMElement* element
          = xml_dom_util::getChildElement(parent, element_name.c_str(), i);
      AminoAcidPtr ptr = std::make_shared<AminoAcid>(element);
      amino_acid_ptr_vec_.push_back(ptr);
      // check if it is an empty acid
      if (ptr->getMonoMass() == 0.0) {
        empty_amino_acid_ptr_ = ptr;
      }
    }
  }
}

AminoAcidPtr AminoAcidBase::getAminoAcidPtrByName(const std::string &name) {
  for (size_t i = 0; i < amino_acid_ptr_vec_.size(); i++) {
    std::string n = amino_acid_ptr_vec_[i]->getName();
    if (n == name) {
      return amino_acid_ptr_vec_[i];
    }
  }
  LOG_DEBUG("Amino acid not found: " + name);
  return AminoAcidPtr(nullptr);
}

AminoAcidPtr AminoAcidBase::getAminoAcidPtrByOneLetter(const std::string &one_letter) {
  for (size_t i = 0; i < amino_acid_ptr_vec_.size(); i++) {
    std::string l = amino_acid_ptr_vec_[i]->getOneLetter();
    if (l == one_letter)  {
      return amino_acid_ptr_vec_[i];
    }
  }
  LOG_DEBUG("Amino acid not found " + one_letter);
  return AminoAcidPtr(nullptr);
}

AminoAcidPtr AminoAcidBase::getAminoAcidPtrByThreeLetter(const std::string &three_letter) {
  for (size_t i = 0; i < amino_acid_ptr_vec_.size(); i++) {
    std::string l = amino_acid_ptr_vec_[i]->getThreeLetter();
    if (l == three_letter) {
      return amino_acid_ptr_vec_[i];
    }
  }
  LOG_DEBUG("Amino acid not found " + three_letter);
  return AminoAcidPtr(nullptr);
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

AminoAcidPtr AminoAcidBase::getAminoAcidPtrFromXml(xercesc::DOMElement * element) {
  std::string name = AminoAcid::getNameFromXml(element);
  AminoAcidPtr acid_ptr = getAminoAcidPtrByName(name);
  return acid_ptr;
}

}  // namespace prot
