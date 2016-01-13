
#include "base/support_peak_type_base.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

SPTypePtrVec SPTypeBase::sp_type_ptr_vec_;
SPTypePtr SPTypeBase::sp_type_ptr_N_TERM_;

void SPTypeBase::initBase(const std::string &file_name){
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    std::string element_name = SupportPeakType::getXmlElementName();
    int prm_peak_type_num = XmlDomUtil::getChildCount(parent, element_name.c_str());
    for (int i = 0; i < prm_peak_type_num; i++) {
      xercesc::DOMElement* element 
          = XmlDomUtil::getChildElement(parent, element_name.c_str(), i);
      SPTypePtr sp_type_ptr(new SupportPeakType(element));
      sp_type_ptr_vec_.push_back(sp_type_ptr);
      if (sp_type_ptr->getName() == getName_N_TERM()) {
        sp_type_ptr_N_TERM_ = sp_type_ptr;
      }
    }
  }
}

SPTypePtr SPTypeBase::getSPTypePtrByName(const std::string &name){
  for (size_t i = 0; i < sp_type_ptr_vec_.size(); i++) {
    std::string n = sp_type_ptr_vec_[i]->getName();
    if (n == name) {
      return sp_type_ptr_vec_[i];
    }
  }
  return SPTypePtr(nullptr);
}

SPTypePtr SPTypeBase::getSPTypePtrById(int id){
  for (size_t i = 0; i < sp_type_ptr_vec_.size(); i++) {
    int n = sp_type_ptr_vec_[i]->getId();
    if (n == id) {
      return sp_type_ptr_vec_[i];
    }
  }
  return SPTypePtr(nullptr);
}

} /* namespace prot */
