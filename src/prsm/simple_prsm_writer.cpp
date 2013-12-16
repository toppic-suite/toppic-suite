/*
 * simple_prsm_writer.cpp
 *
 *  Created on: Dec 9, 2013
 *      Author: xunlikun
 */

#include <iostream>
#include <prsm/simple_prsm_writer.hpp>
#include "base/xml_dom_document.hpp"
#include "base/xml_dom.hpp"

namespace prot {

int SimplePrSMWriter::write(const char * spectrum_file){

		xercesc::DOMImplementation* implementation =  xercesc::DOMImplementationRegistry::getDOMImplementation(X("Core"));
	    xercesc::DOMDocument* doc = implementation->createDocument(0,X("simple_prsm_list"),0);
	    XmlDOMDocument* xml = new XmlDOMDocument(doc);
	    xercesc::DOMElement* root = xml->getDocumentElement();
	    for(unsigned int i = 0;i<matches_.size();i++){
	    	xml->addElement(root,matches_[i]->toXml(xml));
	    }
	    xml->writeXmlDOMDocument(spectrum_file);
//	    delete doc;
//	    delete xml;

	 return 0;
}

void SimplePrSMWriter::addSimplePrSM(SimplePrSMPtrVec matches){
	for(unsigned int i = 0;i<matches.size();i++){
		matches_.push_back(matches[i]);
	}
}

} /* namespace prot */
