#include "xml_dom.hpp"
#include "xml_dom_err_handler.hpp"
 
namespace prot {

static XmlDOMParser* dom_parser;

XmlDOMParser* getXmlDOMInstance() {
    if (dom_parser == NULL) {
        dom_parser = new XmlDOMParser();
    }
    return dom_parser;
}

XmlDOMParser::XmlDOMParser() : parser_(NULL), err_handler_(NULL) {
    xercesc::XMLPlatformUtils::Initialize();
    parser_ = new xercesc::XercesDOMParser();
    err_handler_ = (xercesc::ErrorHandler*) new XmlDOMErrorHandler();
    parser_->setErrorHandler(err_handler_);
}

XmlDOMParser::~XmlDOMParser() {
    if (parser_) {
        delete parser_;
        xercesc::XMLPlatformUtils::Terminate();
        dom_parser = NULL;
    }
}

xercesc::DOMDocument* XmlDOMParser::parse(const char* xml_file) {
  
    parser_->parse(xml_file);
    return parser_->adoptDocument();
}
}
