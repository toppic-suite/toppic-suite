/*
 * prsm_writer.cpp
 *
 *  Created on: Dec 30, 2013
 *      Author: xunlikun
 */

#include "prsm/prsm_writer.hpp"

namespace prot {

PrSMWriter::PrSMWriter(std::string file_name) {
  file_.open(file_name.c_str());
  file_ << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
  file_ << "<prsm_list>" << std::endl;
  XmlDOMImpl* impl = XmlDOMImplFactory::getXmlDOMImplInstance();
  doc_ = new XmlDOMDocument(impl->createDoc("prsm_list"));
  serializer_ = impl->createSerializer();
}

PrSMWriter::~PrSMWriter() {
  file_ << "</prsm_list>" << std::endl;
  file_.close();
  serializer_->release();
  delete doc_;
}

void PrSMWriter::write(PrSM &prsm) {
  xercesc::DOMElement* element = prsm.toXmlElement(doc_);
  std::string str = writeToString(serializer_, element);
  file_ << str << std::endl;
  element->release();
}

} /* namespace prot */
