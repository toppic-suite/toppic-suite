/*
 * prsm_writer.cpp
 *
 *  Created on: Dec 30, 2013
 *      Author: xunlikun
 */

#include <prsm/prsm_writer.hpp>
#include "base/xml_dom_document.hpp"
#include "base/xml_dom.hpp"

namespace prot {
int PrSMWriter::write(const char *prm_file_name){
	xercesc::DOMImplementation* implementation =  xercesc::DOMImplementationRegistry::getDOMImplementation(X("Core"));
	XmlDOMDocument* xml (implementation,"prsm_list");
	xercesc::DOMElement* root = xml->getDocumentElement();
	for(unsigned int i = 0;i<prsms_.size();i++){
		prsms_[i]->appendXml(xml,root);
	}
	xml->writeXmlDOMDocument(prm_file_name);
	delete xml;
}

void PrSMWriter::addSimplePrSM(PrSMPtr matche){
	prsms_.push_back(matche);
}
} /* namespace prot */
