#ifndef XML_WRITER_HPP_
#define XML_WRITER_HPP_

#include <iostream>
#include <fstream>

#include "base/xml_dom_document.hpp"
#include "base/logger.hpp"

namespace prot {

class XmlWriter {
 public:
  XmlWriter(const std::string &file_name,
            const std::string &root);
  ~XmlWriter();
  XmlDOMDocument* getDoc(){return doc_;}
  void write(xercesc::DOMElement* element);
  void close();

 private:
  xercesc::DOMLSSerializer* serializer_;
  XmlDOMDocument* doc_;
  std::ofstream file_;
  std::string root_="";
};

} /* namespace prot */

#endif /* XML_WRITER_HPP_ */
