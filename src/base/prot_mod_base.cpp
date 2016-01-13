#include "base/logger.hpp"
#include "base/prot_mod_base.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

ProtModPtrVec ProtModBase::prot_mod_ptr_vec_;
ProtModPtr ProtModBase::prot_mod_ptr_NONE_;
//ProtModPtr ProtModBase::prot_mod_ptr_NME_;
//ProtModPtr ProtModBase::prot_mod_ptr_NME_ACETYLATION_;

void ProtModBase::initBase(const std::string &file_name) {
  prot::XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    prot::XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    std::string element_name = ProtMod::getXmlElementName();
    int mod_num = XmlDomUtil::getChildCount(parent, element_name.c_str());
    for (int i = 0; i < mod_num; i++) {
      xercesc::DOMElement* element = XmlDomUtil::getChildElement(parent, element_name.c_str(), i);
      ProtModPtr prot_mod_ptr(new ProtMod(element));
      //LOG_DEBUG("index " << i << " nullptr " << (prot_mod_ptr == nullptr));
      prot_mod_ptr_vec_.push_back(prot_mod_ptr);
      if (prot_mod_ptr->getName() == getName_NONE()) {
        prot_mod_ptr_NONE_ = prot_mod_ptr;
      }
      /*
      if (prot_mod_ptr->getName() == getName_NME()) {
        prot_mod_ptr_NME_ = prot_mod_ptr;
      }
      if (prot_mod_ptr->getName() == getName_NME_ACETYLATION()) {
        prot_mod_ptr_NME_ACETYLATION_ = prot_mod_ptr;
      }
      */
    }
  }
}

ProtModPtr ProtModBase::getProtModPtrByName(const std::string &name) {
  for (size_t i = 0; i < prot_mod_ptr_vec_.size(); i++) {
    std::string n = prot_mod_ptr_vec_[i]->getName();
    if (n == name) {
      return prot_mod_ptr_vec_[i];
    }
  }
  LOG_WARN("prot mod nullptr");
  return ProtModPtr(nullptr);
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

ProtModPtr ProtModBase::getProtModPtrFromXml(xercesc::DOMElement * element) {
  std::string name = ProtMod::getNameFromXml(element);
  ProtModPtr prot_mod_ptr = getProtModPtrByName(name);
  return prot_mod_ptr;
}

}
