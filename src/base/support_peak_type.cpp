/*
 * support_peak_type.cpp
 *
 *  Created on: Dec 4, 2013
 *      Author: xunlikun
 */

#include "base/support_peak_type.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

SPTypePtrVec SPTypeFactory::sp_type_ptr_vec_;

void SPTypeFactory::initFactory(const std::string &file_name){
  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMParserInstance();
  if (parser) {
    XmlDOMDocument* doc = new XmlDOMDocument(parser, file_name.c_str());
    if (doc) {
      xercesc::DOMElement* root = doc->getDocumentElement();
      int prm_peak_type_num = getChildCount(root, "support_peak_type");
      for (int i = 0; i < prm_peak_type_num; i++) {
        xercesc::DOMElement* element 
            = getChildElement(root, "support_peak_type", i);
        int id = getIntChildValue(element, "id", 0);
        std::string name = getChildValue(element, "name", 0);
        sp_type_ptr_vec_.push_back(SPTypePtr(new SupportPeakType(id,name)));
      }
      delete doc;
    }
  }
}

SPTypePtr SPTypeFactory::getBaseSPTypePtrByName(const std::string &name){
	for (unsigned int i = 0; i < sp_type_ptr_vec_.size(); i++) {
	    std::string n = sp_type_ptr_vec_[i]->getName();
	    if (n == name) {
	      return sp_type_ptr_vec_[i];
	    }
	  }
	return SPTypePtr(nullptr);
}

SPTypePtr SPTypeFactory::getBaseSPTypePtrById(const int id){
	for (unsigned int i = 0; i < sp_type_ptr_vec_.size(); i++) {
	    int n = sp_type_ptr_vec_[i]->getId();
	    if (n == id) {
	      return sp_type_ptr_vec_[i];
	    }
	}
	return SPTypePtr(nullptr);
}

} /* namespace prot */
