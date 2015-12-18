#include "base/activation.hpp"
#include "base/ion_type_base.hpp"
#include "base/xml_dom_document.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

Activation::Activation(const std::string &name, 
                       IonTypePtr n_ion_type_ptr, 
                       IonTypePtr c_ion_type_ptr):
    name_(name), 
    n_ion_type_ptr_(n_ion_type_ptr),
    c_ion_type_ptr_(c_ion_type_ptr) {
    }

Activation::Activation(xercesc::DOMElement * element) {
  name_ = XmlDomUtil::getChildValue(element, "name", 0);
  std::string ion_type_name = XmlDomUtil::getChildValue(element, "n_ion_type",0);
  n_ion_type_ptr_ = IonTypeBase::getIonTypePtrByName(ion_type_name);
  ion_type_name = XmlDomUtil::getChildValue(element, "c_ion_type", 0);
  c_ion_type_ptr_ = IonTypeBase::getIonTypePtrByName(ion_type_name);
}

void Activation::appendNameToXml(XmlDOMDocument* xml_doc,
                                 xercesc::DOMElement* parent){
  std::string element_name = Activation::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  xml_doc->addElement(element, "name", name_.c_str());
  parent->appendChild(element);
}

std::string Activation::getNameFromXml(xercesc::DOMElement * element) {
  std::string name = XmlDomUtil::getChildValue(element, "name", 0);
  return name;
}

} 
