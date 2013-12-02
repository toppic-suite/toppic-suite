#ifndef PROT_XML_DOM_DOCUMENT_HPP_
#define PROT_XML_DOM_DOCUMENT_HPP_

#include <string>

#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/sax/HandlerBase.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/PlatformUtils.hpp>

#include "base/xml_dom.hpp"

namespace prot {

class XmlDOMDocument {

 public:
  XmlDOMDocument(XmlDOMParser* parser, const char* xml_file);
  ~XmlDOMDocument();

  xercesc::DOMElement* createElement(const char* tag);

  xercesc::DOMText* createTextNode(const char* text);

  void addElement(xercesc::DOMElement* element, 
                  const char* tag, const char* value);

  xercesc::DOMElement* getDocumentElement() {
    return doc_->getDocumentElement();
  }

 private:
  xercesc::DOMDocument* doc_;
};

xercesc::DOMElement* getChildElement(xercesc::DOMElement* parent, 
                                const char* tag, int index);

std::string getChildValue(xercesc::DOMElement* parent, 
                          const char* child_tag, int index);

double getDoubleChildValue(xercesc::DOMElement* parent, 
                          const char* child_tag, int index);

int getIntChildValue(xercesc::DOMElement* parent,  
                     const char* child_tag, int index);

bool getBoolChildValue(xercesc::DOMElement* parent,
                       const char* child_tag, int index);

int getChildCount(xercesc::DOMElement * parent, const char* child_tag);

std::string getAttributeValue(xercesc::DOMElement * parent,
                                const char* attribute_tag);

std::string convertToString(double value);

std::string convertToString(bool value);

}
#endif
