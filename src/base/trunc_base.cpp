#include "base/logger.hpp"
#include "base/trunc_base.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

TruncPtrVec TruncBase::trunc_ptr_vec_;

void TruncBase::initBase(const std::string &file_name) {
  prot::XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    prot::XmlDOMDocument doc(parser, file_name.c_str());
    xercesc::DOMElement* parent = doc.getDocumentElement();
    std::string element_name = Trunc::getXmlElementName();
    int trunc_num = XmlDomUtil::getChildCount(parent, element_name.c_str());
    for (int i = 0; i < trunc_num; i++) {
      xercesc::DOMElement* element 
          = XmlDomUtil::getChildElement(parent, element_name.c_str(), i);
      TruncPtr trunc_ptr(new Trunc(element));
      trunc_ptr_vec_.push_back(trunc_ptr);
    }
  }
}

TruncPtr TruncBase::getTruncPtrByName(const std::string &name) {
  for (size_t i = 0; i < trunc_ptr_vec_.size(); i++) {
    std::string n = trunc_ptr_vec_[i]->getName();
    if (n == name) {
      return trunc_ptr_vec_[i];
    }
  }
  return TruncPtr(nullptr);
}

TruncPtr TruncBase::getTruncPtrFromXml(xercesc::DOMElement * element) {
  std::string name = Trunc::getNameFromXml(element);
  TruncPtr trunc_ptr = getTruncPtrByName(name);
  return trunc_ptr;
}

}
