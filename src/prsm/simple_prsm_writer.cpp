/*
 * simple_prsm_writer.cpp
 *
 *  Created on: Dec 9, 2013
 *      Author: xunlikun
 */

#include <iostream>
#include "prsm/simple_prsm_writer.hpp"
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
    for(unsigned int i=0;i<simple_prsms.size();i++){
      xercesc::DOMElement* element = simple_prsms[i]->toXml(doc_);
      std::string str = writeToString(serializer_, element);
      writeToStreamByRemovingDoubleLF(file_, str);
      element->release();
    }
}

} /* namespace prot */
