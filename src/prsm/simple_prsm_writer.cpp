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
SimplePrSMWriter::SimplePrSMWriter(std::string file_name){
	file_.open(file_name.c_str());
	file_ << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
	file_ << "<simple_prsm_list>" << std::endl;
	XmlDOMImpl* impl = XmlDOMImplFactory::getXmlDOMImplInstance();
	doc_ = new XmlDOMDocument(impl->createDoc("simple_prsm_list"));
	serializer_ = impl->createSerializer();
}
SimplePrSMWriter::~SimplePrSMWriter(){
	file_ << "</simple_prsm_list>" << std::endl;
	file_.close();
	serializer_->release();
	delete doc_;
}
void SimplePrSMWriter::write(SimplePrSMPtrVec simple_prsms){
	for(int i=0;i<simple_prsms.size();i++){
	  xercesc::DOMElement* element = simple_prsms[i]->toXml(doc_);//>toXmlElement(doc_);
	  std::string str = writeToString(serializer_, element);
	  writeToStreamByRemovingDoubleLF(file_, str);
	  //file_ << str << std::endl;
	  element->release();
	}
}

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
