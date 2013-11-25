#ifndef PROT_XML_DOM_DOCUMENT_H_
#define PROT_XML_DOM_DOCUMENT_H_

#include <string>

#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/sax/HandlerBase.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/PlatformUtils.hpp>

#include "xml_dom.hpp"

namespace prot {

class XmlDOMDocument {
       
public:
    XmlDOMDocument(XmlDOMParser* parser, const char* xml_file);
    ~XmlDOMDocument();

    xercesc::DOMElement* getElement(const char* tag, int index);

    std::string getChildValue(const char* parent_tag, int parent_index, 
            const char* child_tag);
    
    std::string getAttributeValue(const char* element_tag,  
            int element_index, 
            const char* attribute_tag);

    int getChildCount(const char* parent_tag, int parent_index, 
            const char* child_tag);

private:
    xercesc::DOMDocument* doc_;
    XmlDOMDocument();
    XmlDOMDocument(const XmlDOMDocument&); 
};

std::string getChildValue(xercesc::DOMElement* parent,  
        const char* child_tag);

double getDoubleChildValue(xercesc::DOMElement* parent,  
        const char* child_tag);

int getIntChildValue(xercesc::DOMElement* parent,  
        const char* child_tag);
}
#endif
