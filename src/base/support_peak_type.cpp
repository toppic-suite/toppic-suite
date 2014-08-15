
#include "base/support_peak_type.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

SPTypePtrVec SPTypeFactory::sp_type_ptr_vec_;

SupportPeakType::SupportPeakType(int id, const std::string &name) {
  id_=id;
  name_=name;
}

void SPTypeFactory::initFactory(const std::string &file_name){
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    int prm_peak_type_num = getChildCount(parent, "support_peak_type");
    for (int i = 0; i < prm_peak_type_num; i++) {
      xercesc::DOMElement* element 
          = getChildElement(parent, "support_peak_type", i);
      int id = getIntChildValue(element, "id", 0);
      std::string name = getChildValue(element, "name", 0);
      sp_type_ptr_vec_.push_back(SPTypePtr(new SupportPeakType(id,name)));
    }
  }
}

SPTypePtr SPTypeFactory::getBaseSPTypePtrByName(const std::string &name){
  for (size_t i = 0; i < sp_type_ptr_vec_.size(); i++) {
    std::string n = sp_type_ptr_vec_[i]->getName();
    if (n == name) {
      return sp_type_ptr_vec_[i];
    }
  }
  return SPTypePtr(nullptr);
}

SPTypePtr SPTypeFactory::getBaseSPTypePtrById(int id){
  for (size_t i = 0; i < sp_type_ptr_vec_.size(); i++) {
    int n = sp_type_ptr_vec_[i]->getId();
    if (n == id) {
      return sp_type_ptr_vec_[i];
    }
  }
  return SPTypePtr(nullptr);
}

} /* namespace prot */
