#include "base/xml_dom_err_handler.hpp"

namespace prot {

void XmlDOMErrorHandler::fatalError(const xercesc::SAXParseException &ex) {
  std::cerr << "Fatal parsing error at line" << ex.getLineNumber() << "\n"; 
  exit(-1);
}

}
