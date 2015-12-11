#include "base/logger.hpp"
#include "base/ptm_base.hpp"
#include "base/trunc_base.hpp"
#include "base/mod_base.hpp"
#include "base/prot_mod.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

ProtMod::ProtMod(const std::string &name, const std::string &type,
                 TruncPtr trunc_ptr, ModPtr mod_ptr): 
    name_(name),
    type_(type),
    trunc_ptr_(trunc_ptr),
    mod_ptr_(mod_ptr) {
      mod_pos_ = trunc_ptr->getTruncLen();
      prot_shift_ = trunc_ptr_->getShift() + mod_ptr_->getShift();
      pep_shift_ = mod_ptr_->getShift();
    }

ProtMod::ProtMod(xercesc::DOMElement* element) { 
  name_ = XmlDomUtil::getChildValue(element, "name", 0);
  type_ = XmlDomUtil::getChildValue(element, "type", 0);
  std::string trunc_element_name = Trunc::getXmlElementName();
  xercesc::DOMElement* trunc_element 
      = XmlDomUtil::getChildElement(element, trunc_element_name.c_str(), 0);
  trunc_ptr_ = TruncBase::getTruncPtrFromXml(trunc_element);
  std::string mod_element_name = Mod::getXmlElementName();
  xercesc::DOMElement* mod_element 
      = XmlDomUtil::getChildElement(element, mod_element_name.c_str(), 0);
  mod_ptr_= ModBase::getModPtrFromXml(mod_element); 
  mod_pos_ = trunc_ptr_->getTruncLen();
  prot_shift_ = trunc_ptr_->getShift() + mod_ptr_->getShift();
  pep_shift_ = mod_ptr_->getShift();
}

void ProtMod::appendNameToXml(XmlDOMDocument* xml_doc,
                              xercesc::DOMElement* parent){
  std::string element_name = ProtMod::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  xml_doc->addElement(element, "name", name_.c_str());
  parent->appendChild(element);
}

std::string ProtMod::getNameFromXml(xercesc::DOMElement * element) {
  std::string name = XmlDomUtil::getChildValue(element, "name", 0);
  return name;
}

}
