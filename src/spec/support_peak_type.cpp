/*
 * support_peak_type.cpp
 *
 *  Created on: Dec 4, 2013
 *      Author: xunlikun
 */

#include "spec/support_peak_type.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

SupportPeakTypePtrVec getSupportPeakTypePtrVecInstance(const char* file_name){
	SupportPeakTypePtrVec support_peak_type_list;
		  XmlDOMParser* parser = XmlDOMParserFactory::getXmlDOMInstance();
		  if (parser) {
		    XmlDOMDocument* doc = new XmlDOMDocument(parser, file_name);
		    if (doc) {
		      xercesc::DOMElement* root = doc->getDocumentElement();
		      xercesc::DOMElement* parent = getChildElement(root, "support_peak_type_list", 0);
		      int prm_peak_type_num = getChildCount(parent, "support_peak_type");
		      for (int i = 0; i < prm_peak_type_num; i++) {
		        xercesc::DOMElement* element = getChildElement(parent, "support_peak_type", i);
		        std::string name = getChildValue(element, "name", 0);
		        int id = getIntChildValue(element, "id", 0);
		        support_peak_type_list.push_back(SupportPeakTypePtr(
		                new SupportPeakType(id,name)));

		      }
		      delete doc;
		    }
		    delete parser;
		  }
		  return support_peak_type_list;
}

SupportPeakTypePtr getSupportPeakTypePtrByName(SupportPeakTypePtrVec &support_peak_type_list,
                         const std::string &name){
	for (unsigned int i = 0; i < support_peak_type_list.size(); i++) {
	    std::string n = support_peak_type_list[i]->getName();
	    if (n == name) {
	      return support_peak_type_list[i];
	    }
	  }
	return SupportPeakTypePtr(nullptr);
}

SupportPeakTypePtr getSupportPeakTypePtrById(SupportPeakTypePtrVec &support_peak_type_list,
                         const int id){
	for (unsigned int i = 0; i < support_peak_type_list.size(); i++) {
	    int n = support_peak_type_list[i]->getId();
	    if (n == id) {
	      return support_peak_type_list[i];
	    }
	}
	return SupportPeakTypePtr(nullptr);
}

} /* namespace prot */
