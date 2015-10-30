#include "base/activation_base.hpp"

namespace prot {

ActivationPtrVec ActivationBase::activation_ptr_vec_;

/* Methods for ActivationFactory */
void ActivationBase::initBase(const std::string &file_name){
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

} /* namespace prot */
