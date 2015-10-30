/*
 * author  Xiaowen Liu
 * date    2013-11-1
 */
#include <base/logger.hpp>

#include "base/acid_base.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

AcidPtrVec AcidBase::acid_ptr_vec_; 
AcidPtr AcidBase::empty_acid_ptr_;

void AcidBase::initBase(const std::string &file_name) {
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    int acid_num = getChildCount(parent, "amino_acid");
    LOG_DEBUG("acid num " << acid_num);
    for (int i = 0; i < acid_num; i++) {
      xercesc::DOMElement* element = getChildElement(parent, "amino_acid", i);
      std::string name = getChildValue(element, "name", 0);
      std::string one_letter = getChildValue(element, "one_letter", 0);
      std::string three_letter = getChildValue(element, "three_letter", 0);
      std::string composition = getChildValue(element, "composition", 0);
      double mono_mass = getDoubleChildValue(element, "mono_mass", 0);
      double avg_mass = getDoubleChildValue(element, "average_mass", 0);
      AcidPtr ptr(new Acid(name, one_letter, three_letter, 
                           composition, mono_mass, avg_mass));
      acid_ptr_vec_.push_back(ptr);
      // check if it is an empty acid
      if (ptr->getMonoMass() == 0.0) {
        empty_acid_ptr_ = ptr;
      }
    }
  }
}

/**
 * Returns an amino acid based on the the name. Returns null if the amino
 * acid name does not exist.
 */
AcidPtr AcidBase::getBaseAcidPtrByName(const std::string &name) {
  for (size_t i = 0; i < acid_ptr_vec_.size(); i++) {
    std::string n = acid_ptr_vec_[i]->getName();
    if (n == name) {
      return acid_ptr_vec_[i];
    }
  }
  return AcidPtr(nullptr);
}

/**
 * Returns an amino acid based on the one letter representation. Returns
 * null if the one letter representation does not exist.
 */
AcidPtr AcidBase::getBaseAcidPtrByOneLetter(const std::string &one_letter) {
  for (size_t i = 0; i < acid_ptr_vec_.size(); i++) {
    std::string l = acid_ptr_vec_[i]->getOneLetter();
    if (l == one_letter)  {
      return acid_ptr_vec_[i];
    }
  }
  LOG_DEBUG( "Acid not found " + one_letter);
  return AcidPtr(nullptr);
}

/**
 * Returns an amino acid based on the three letter representation. Returns
 * null if the three letter representation does not exist.
 */
AcidPtr AcidBase::getBaseAcidPtrByThreeLetter(const std::string &three_letter) {
  for (size_t i = 0; i < acid_ptr_vec_.size(); i++) {
    std::string l = acid_ptr_vec_[i]->getThreeLetter();
    if (l == three_letter) {
      return acid_ptr_vec_[i];
    }
  }
  LOG_DEBUG( "Acid not found " + three_letter);
  return AcidPtr(nullptr);
}

/**
 * Checks if the list contains an amino acid with the specific name.
 */
bool AcidBase::baseContainsName(const std::string &name) {
  return getBaseAcidPtrByName(name).get() != nullptr;
}

/**
 * Checks if the list contains an amino acid with the specific one letter
 * representation.
 */
bool AcidBase::baseContainsOneLetter(const std::string &one_letter) {
  return getBaseAcidPtrByOneLetter(one_letter).get() != nullptr;
}

/**
 * Checks if the list contains an amino acid with the specific three letter
 * representation.
 */
bool AcidBase::baseContainsThreeLetter(const std::string &three_letter) {
  return getBaseAcidPtrByThreeLetter(three_letter).get() != nullptr;
}

/**
 * Converts a protein sequence (with one letter representation of amino
 * acids) to an amino acid array.
 */
AcidPtrVec AcidBase::convertSeqToAcidSeq(const std::string &seq) {
  AcidPtrVec acid_seq;
  if (seq.length() > 0) {
    for (size_t i = 0; i < seq.length(); i++) {
      acid_seq.push_back(getBaseAcidPtrByOneLetter(seq.substr(i, 1)));
    }
  }
  return acid_seq;
}

} /* end namespace */

