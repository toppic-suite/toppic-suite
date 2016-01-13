#ifndef PROT_BASE_XML_DOM_ERR_HANDLER_HPP_
#define PROT_BASE_XML_DOM_ERR_HANDLER_HPP_

#include <iostream>

#include <xercesc/sax/HandlerBase.hpp>

namespace prot {

class XmlDOMErrorHandler : public xercesc::HandlerBase {
 public:
  void fatalError(const xercesc::SAXParseException &ex);
};
 
}
#endif
