#ifndef PROTOMICS_XML_DOM_ERR_HANDLER_H_
#define PROTOMICS_XML_DOM_ERR_HANDLER_H_

#include <iostream>

#include <xercesc/sax/HandlerBase.hpp>

namespace proteomics {

class XmlDOMErrorHandler : public xercesc::HandlerBase {
    public:
    void fatalError(const xercesc::SAXParseException &ex) {
        std::cerr << "Fatal parsing error at line" << (int)ex.getLineNumber() << "\n"; 
        exit(-1);
    }
};
 
}
#endif
