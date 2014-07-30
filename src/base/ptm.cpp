/*
 * author  Xiaowen Liu
 * date    2013-11-1
 */
#include <iostream>

#include "base/ptm.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

PtmPtrVec PtmFactory::ptm_ptr_vec_;

Ptm::Ptm(const std::string &abbr_name, 
         double mono_mass) {
  abbr_name_ = abbr_name;
  mono_mass_ = mono_mass;
}

bool Ptm::isEmpty() {
  if (mono_mass_ == 0.0) {
    return true;
  }
  else {
    return false;
  }
}

bool Ptm::isAcetylation() {
  if (abbr_name_ == PTM_ACETYLATION) {
    return true;
  }
  else {
    return false;
  }
}

void Ptm::appendxml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
  xercesc::DOMElement* element = xml_doc->createElement("modification");
  xml_doc->addElement(element, "abbr_name", abbr_name_.c_str());
  std::string str = convertToString(mono_mass_);
  xml_doc->addElement(element, "mono_mass", str.c_str());
  parent->appendChild(element);
}


void PtmFactory::initFactory(const std::string &file_name) {
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    int ptm_num = getChildCount(parent, "modification");
    for (int i = 0; i < ptm_num; i++) {
      xercesc::DOMElement* element = getChildElement(parent, "modification", i);
      std::string abbr_name = getChildValue(element, "abbreviation", 0);
      double mono_mass = getDoubleChildValue(element, "mono_mass", 0);
      ptm_ptr_vec_.push_back(PtmPtr(new Ptm(abbr_name, mono_mass)));

    }
  }
}

PtmPtr PtmFactory::findEmptyPtmPtr() {
  for (size_t i = 0; i < ptm_ptr_vec_.size(); i++) {
    if (ptm_ptr_vec_[i]->isEmpty()) {
      return ptm_ptr_vec_[i];
    }
  }
  throw "Empty ptm does not exist!";
}

/**
 *   Returns a PTM based on the abbreviation name. Returns null if the
 *   abbreviation name does not exist.
 */
PtmPtr PtmFactory::getBasePtmPtrByAbbrName(const std::string &abbr_name) {
  for (size_t i = 0; i < ptm_ptr_vec_.size(); i++) {
    std::string n = ptm_ptr_vec_[i]->getAbbrName();
    if (n == abbr_name) {
      return ptm_ptr_vec_[i];
    }
  }
  return PtmPtr(nullptr);
}

/**
 * Checks if the list contains an amino acid with the specific name.
 */
bool PtmFactory::baseContainAbbrName(const std::string &abbr_name) {
  return getBasePtmPtrByAbbrName(abbr_name).get() != nullptr;
}

PtmPtr PtmFactory::addBasePtm(const std::string &abbr_name, double mono_mass) {
  PtmPtr ptm_ptr = getBasePtmPtrByAbbrName(abbr_name);
  if (ptm_ptr.get() == nullptr) {
    PtmPtr new_ptm(new Ptm(abbr_name, mono_mass));
    ptm_ptr_vec_.push_back(new_ptm);
    return new_ptm;
  }
  else {
    return ptm_ptr;
  }
}

}

