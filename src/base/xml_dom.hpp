#ifndef PROT_XML_DOM_HPP_
#define PROT_XML_DOM_HPP_

#include <string>

#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/util/PlatformUtils.hpp>

namespace prot {

class XmlDOMParser {
public:
    XmlDOMParser();
    ~XmlDOMParser();

    xercesc::DOMDocument* parse(const char* xml_file);

private:
    xercesc::XercesDOMParser* parser_;
    xercesc::ErrorHandler*    err_handler_;
};

XmlDOMParser* getXmlDOMInstance();
 
}
#endif
