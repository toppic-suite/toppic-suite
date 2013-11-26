#ifndef PROT_XML_DOM_ERR_HANDLER_HPP_
#define PROT_XML_DOM_ERR_HANDLER_HPP_

#include <iostream>

#include <xercesc/sax/HandlerBase.hpp>

namespace prot {

class XmlDOMErrorHandler : public xercesc::HandlerBase {
    public:
    void fatalError(const xercesc::SAXParseException &ex) {
        std::cerr << "Fatal parsing error at line" << (int)ex.getLineNumber() << "\n"; 
        exit(-1);
    }
};
 
}
#endif
