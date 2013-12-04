/*
 * prm_peak_type.cpp
 *
 *  Created on: Dec 4, 2013
 *      Author: xunlikun
 */

#include "spec/prm_peak_type.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

PrmPeakTypePtrVec getPrmPeakTypePtrVecInstance(const char* file_name){
	PrmPeakTypePtrVec prm_peak_type_list;
	  XmlDOMParser* parser = getXmlDOMInstance();
	  if (parser) {
	    XmlDOMDocument* doc = new XmlDOMDocument(parser, file_name);
	    if (doc) {
	      xercesc::DOMElement* root = doc->getDocumentElement();
	      xercesc::DOMElement* parent = getChildElement(root, "prm_peak_type_list", 0);
	      int prm_peak_type_num = getChildCount(parent, "prm_peak_type");
	      for (int i = 0; i < prm_peak_type_num; i++) {
	        xercesc::DOMElement* element = getChildElement(parent, "prm_peak_type", i);
	        std::string name = getChildValue(element, "name", 0);
	        int id = getIntChildValue(element, "id", 0);
	        prm_peak_type_list.push_back(PrmPeakTypePtr(
	                new PrmPeakType(id,name)));

	      }
	      delete doc;
	    }
	    delete parser;
	  }
	  return prm_peak_type_list;
}

PrmPeakTypePtr getPrmPeakTypePtrByName(PrmPeakTypePtrVec &prm_peak_type_list,
                         const std::string &name){
	for (unsigned int i = 0; i < prm_peak_type_list.size(); i++) {
	    std::string n = prm_peak_type_list[i]->getName();
	    if (n == name) {
	      return prm_peak_type_list[i];
	    }
	  }
	  return PrmPeakTypePtr(nullptr);
}

PrmPeakTypePtr getPrmPeakTypePtrById(PrmPeakTypePtrVec &prm_peak_type_list,
                         const int id ){
	for (unsigned int i = 0; i < prm_peak_type_list.size(); i++) {
	    int n = prm_peak_type_list[i]->getId();
	    if (n == id) {
	      return prm_peak_type_list[i];
	    }
	  }
	  return PrmPeakTypePtr(nullptr);
}

} /* namespace prot */
