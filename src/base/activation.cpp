
#include "base/activation.hpp"
#include "base/ion_type.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

ActivationPtrVec ActivationFactory::activation_ptr_vec_;

Activation::Activation(const std::string &name, 
                       IonTypePtr n_ion_type_ptr, 
                       IonTypePtr c_ion_type_ptr) {
  name_ = name;
  n_ion_type_ptr_ = n_ion_type_ptr;
  c_ion_type_ptr_ = c_ion_type_ptr;
}

Activation::Activation(xercesc::DOMElement * element) {
  name_ = getChildValue(element, "name", 0);
  std::string ion_type_name = getChildValue(element, "n_ion_type",0);
  n_ion_type_ptr_ = IonTypeFactory::getBaseIonTypePtrByName(ion_type_name);
  ion_type_name = getChildValue(element, "c_ion_type", 0);
  c_ion_type_ptr_ = IonTypeFactory::getBaseIonTypePtrByName(ion_type_name);
}

void Activation::appendXml(XmlDOMDocument* xml_doc,
                           xercesc::DOMElement* parent){
    xercesc::DOMElement* element = xml_doc->createElement("activation");
    xml_doc->addElement(element, "name", name_.c_str());
    parent->appendChild(element);
}

/* Methods for ActivationFactory */
void ActivationFactory::initFactory(const std::string &file_name){
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    int activation_num = getChildCount(parent, "activation");
    for (int i = 0; i < activation_num; i++) {
      xercesc::DOMElement* element 
          = getChildElement(parent, "activation", i);
      ActivationPtr ptr(new Activation(element));
      activation_ptr_vec_.push_back(ptr);
    }
  }
}

ActivationPtr ActivationFactory::getBaseActivationPtrByName(
    const std::string &name){
  for (size_t i = 0; i < activation_ptr_vec_.size(); i++) {
    std::string n = activation_ptr_vec_[i]->getName();
    if (n == name) {
      return activation_ptr_vec_[i];
    }
  }
  return ActivationPtr(nullptr);
}

} /* namespace prot */
