#include "base/logger.hpp"
#include "base/mod.hpp"
#include "base/residue_base.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

Mod::Mod(ResiduePtr ori_residue_ptr, ResiduePtr mod_residue_ptr):
    ori_residue_ptr_(ori_residue_ptr),
    mod_residue_ptr_(mod_residue_ptr) {
    }

Mod::Mod(xercesc::DOMElement* element) { 
  xercesc::DOMElement* ori_residue_element 
      = XmlDomUtil::getChildElement(element, "ori_residue", 0);
  ResiduePtr tmp_ori_ptr(new Residue(ori_residue_element)); 
  ori_residue_ptr_ = ResidueBase::getBaseResiduePtr(tmp_ori_ptr);
  xercesc::DOMElement* mod_residue_element 
      = XmlDomUtil::getChildElement(element, "mod_residue", 0);
  ResiduePtr tmp_mod_ptr(new Residue(mod_residue_element)); 
  mod_residue_ptr_ = ResidueBase::getBaseResiduePtr(tmp_mod_ptr);
}

void Mod::appendToXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent){
  std::string element_name = Mod::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  ori_residue_ptr_->appendXml(xml_doc, element, "ori_residue");
  mod_residue_ptr_->appendXml(xml_doc, element, "mod_residue");
  parent->appendChild(element);
}

}
