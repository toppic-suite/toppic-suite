/*
 * author  Xiaowen Liu
 * date    2013-11-1
 */
#include "base/logger.hpp"
#include "base/acid_base.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

AcidPtrVec AcidBase::acid_ptr_vec_; 
AcidPtr AcidBase::empty_acid_ptr_;

void AcidBase::initBase(const std::string &file_name) {
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    std::string element_name = Acid::getXmlElementName();
    int acid_num = XmlDomUtil::getChildCount(parent, element_name.c_str());
    LOG_DEBUG("acid num " << acid_num);
    for (int i = 0; i < acid_num; i++) {
      xercesc::DOMElement* element 
          = XmlDomUtil::getChildElement(parent, element_name.c_str(), i);
      AcidPtr ptr(new Acid(element));
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
AcidPtr AcidBase::getAcidPtrByName(const std::string &name) {
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
AcidPtr AcidBase::getAcidPtrByOneLetter(const std::string &one_letter) {
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
AcidPtr AcidBase::getAcidPtrByThreeLetter(const std::string &three_letter) {
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
bool AcidBase::containsName(const std::string &name) {
  return getAcidPtrByName(name).get() != nullptr;
}

/**
 * Checks if the list contains an amino acid with the specific one letter
 * representation.
 */
bool AcidBase::containsOneLetter(const std::string &one_letter) {
  return getAcidPtrByOneLetter(one_letter).get() != nullptr;
}

/**
 * Checks if the list contains an amino acid with the specific three letter
 * representation.
 */
bool AcidBase::containsThreeLetter(const std::string &three_letter) {
  return getAcidPtrByThreeLetter(three_letter).get() != nullptr;
}

AcidPtr AcidBase::getAcidPtrFromXml(xercesc::DOMElement * element) {
  std::string name = Acid::getNameFromXml(element);
  AcidPtr acid_ptr = getAcidPtrByName(name);
  return acid_ptr;
}

} /* end namespace */

