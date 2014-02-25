/*
 * xml_writer.hpp
 *
 *  Created on: Feb 24, 2014
 *      Author: xunlikun
 */

#ifndef XML_WRITER_HPP_
#define XML_WRITER_HPP_

#include <iostream>
#include <fstream>

#include "base/xml_dom_document.hpp"
#include "base/logger.hpp"

namespace prot {

class XmlWriter {
 public:
  XmlWriter(std::string file_name,std::string root);
  ~XmlWriter();
  XmlDOMDocument* getDoc(){return doc_;}
  void write(xercesc::DOMElement* element);
 private:
  xercesc::DOMLSSerializer* serializer_;
  XmlDOMDocument* doc_;
  std::ofstream file_;
  std::string root_="";
};

} /* namespace prot */

#endif /* XML_WRITER_HPP_ */
