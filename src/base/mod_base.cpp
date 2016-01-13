#include <iostream>
#include <vector>
#include <memory>
#include <algorithm>

#include "base/ptm_base.hpp"
#include "base/mod_base.hpp"
#include "base/logger.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

ModPtrVec ModBase::mod_ptr_vec_;
ModPtr ModBase::none_mod_ptr_;
ModPtr ModBase::c57_mod_ptr_;
ModPtr ModBase::c58_mod_ptr_;

void ModBase::initBase(const std::string &file_name) {
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    std::string element_name = Mod::getXmlElementName();
    int mod_num = XmlDomUtil::getChildCount(parent, element_name.c_str());
    for (int i = 0; i < mod_num; i++) {
      xercesc::DOMElement* element = XmlDomUtil::getChildElement(parent, element_name.c_str(), i);
      ModPtr mod_ptr(new Mod(element));
      mod_ptr_vec_.push_back(mod_ptr);
      // check empty ptr
      if (mod_ptr->getOriResiduePtr() == ResidueBase::getEmptyResiduePtr() 
          && mod_ptr->getModResiduePtr() ==ResidueBase::getEmptyResiduePtr()) {
        none_mod_ptr_ = mod_ptr;
      }
      if (mod_ptr->getModResiduePtr()->getAcidPtr()->getOneLetter() == "C"
          && mod_ptr->getModResiduePtr()->getPtmPtr() == PtmBase::getPtmPtr_C57()) {
        c57_mod_ptr_ = mod_ptr;
      }
      if (mod_ptr->getModResiduePtr()->getAcidPtr()->getOneLetter() == "C"
          && mod_ptr->getModResiduePtr()->getPtmPtr() == PtmBase::getPtmPtr_C58()) {
        c58_mod_ptr_ = mod_ptr;
      }
    }
    if (none_mod_ptr_ == nullptr || c57_mod_ptr_ == nullptr || c58_mod_ptr_ == nullptr) {
      LOG_WARN("mod missing!");
    }
  }
}

ModPtr ModBase::getBaseModPtr(ModPtr mod_ptr) {
  for (size_t i = 0; i < mod_ptr_vec_.size(); i++) {
    if (mod_ptr_vec_[i]->isSame(mod_ptr)) {
      return mod_ptr_vec_[i];
    }
  }
  mod_ptr_vec_.push_back(mod_ptr);
  return mod_ptr;
}

ModPtr ModBase::getModPtrFromXml(xercesc::DOMElement * element) {
  ModPtr ptr(new Mod(element));
  return getBaseModPtr(ptr);
}
}

