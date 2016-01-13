#include "base/ion_type_base.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

IonTypePtrVec IonTypeBase::ion_type_ptr_vec_; 
IonTypePtr IonTypeBase::ion_type_ptr_B_; 
IonTypePtr IonTypeBase::ion_type_ptr_PREC_; 

void IonTypeBase::initBase(const std::string &file_name){
  prot::XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    prot::XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    std::string element_name = IonType::getXmlElementName();
    int ion_type_num = XmlDomUtil::getChildCount(parent, element_name.c_str());
    for (int i = 0; i < ion_type_num; i++) {
      xercesc::DOMElement* element = XmlDomUtil::getChildElement(parent, element_name.c_str(), i);
      IonTypePtr ion_type_ptr(new IonType(element));
      ion_type_ptr_vec_.push_back(ion_type_ptr);
      if (ion_type_ptr->getName() == getName_B()) {
        ion_type_ptr_B_ = ion_type_ptr;
      }
      if (ion_type_ptr->getName() == getName_PREC()) {
        ion_type_ptr_PREC_ = ion_type_ptr;
      }
    }
  }
}

IonTypePtr IonTypeBase::getIonTypePtrByName(const std::string &name){
  for (size_t i = 0; i < ion_type_ptr_vec_.size(); i++) {
    std::string n = ion_type_ptr_vec_[i]->getName();
    if (n == name) {
      return ion_type_ptr_vec_[i];
    }
  }
  return IonTypePtr(nullptr);
}

}
