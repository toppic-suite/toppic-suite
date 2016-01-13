#include "base/activation_base.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

ActivationPtrVec ActivationBase::activation_ptr_vec_;

void ActivationBase::initBase(const std::string &file_name){
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    std::string element_name = Activation::getXmlElementName();
    int activation_num = XmlDomUtil::getChildCount(parent, element_name.c_str());
    for (int i = 0; i < activation_num; i++) {
      xercesc::DOMElement* element 
          = XmlDomUtil::getChildElement(parent, element_name.c_str(), i);
      ActivationPtr ptr(new Activation(element));
      activation_ptr_vec_.push_back(ptr);
    }
  }
}

ActivationPtr ActivationBase::getActivationPtrByName(
    const std::string &name){
  for (size_t i = 0; i < activation_ptr_vec_.size(); i++) {
    std::string n = activation_ptr_vec_[i]->getName();
    if (n == name) {
      return activation_ptr_vec_[i];
    }
  }
  return ActivationPtr(nullptr);
}

ActivationPtr ActivationBase::getActivationPtrFromXml(xercesc::DOMElement * element) {
  std::string name = Activation::getNameFromXml(element);
  ActivationPtr activation_ptr = getActivationPtrByName(name);
  return activation_ptr;
}

} 
