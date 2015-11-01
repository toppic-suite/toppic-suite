#include "base/logger.hpp"
#include "base/ptm_base.hpp"
#include "base/trunc_base.hpp"
#include "base/prot_mod.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

ProtMod::ProtMod(const std::string &name, TruncPtr trunc_ptr, 
                 PtmPtr ptm_ptr): 
    name_(name),
    trunc_ptr_(trunc_ptr),
    ptm_ptr_(ptm_ptr) {
      prot_shift_ = trunc_ptr_->getShift() + ptm_ptr_->getMonoMass();
      pep_shift_ = ptm_ptr_->getMonoMass();
    }

ProtMod::ProtMod(xercesc::DOMElement* element) { 
  name_ = XmlDomUtil::getChildValue(element, "name", 0);
  std::string trunc_element_name = Trunc::getXmlElementName();
  xercesc::DOMElement* trunc_element 
      = XmlDomUtil::getChildElement(element, trunc_element_name.c_str(), 0);
  trunc_ptr_ = TruncBase::getTruncPtrFromXml(trunc_element);
  std::string ptm_element_name = Ptm::getXmlElementName();
  xercesc::DOMElement* ptm_element 
      = XmlDomUtil::getChildElement(element, ptm_element_name.c_str(), 0);
  ptm_ptr_ = PtmBase::getPtmPtrFromXml(ptm_element);
  prot_shift_ = trunc_ptr_->getShift() + ptm_ptr_->getMonoMass();
  pep_shift_ = ptm_ptr_->getMonoMass();
}

void ProtMod::appendNameToXml(XmlDOMDocument* xml_doc,
                              xercesc::DOMElement* parent){
  std::string element_name = ProtMod::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  xml_doc->addElement(element, "name", name_.c_str());
  parent->appendChild(element);
}

}
