#include <stdlib.h>
#include <string>
#include <sstream>
#include <exception>
#include <algorithm>

#include "base/logger.hpp"
#include "base/xml_dom.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

XmlDOMDocument::XmlDOMDocument(XmlDOMParser* parser, 
                               const char* xml_file) : doc_(NULL) {
  try {
    doc_ = parser->parse(xml_file);
  }
  catch (std::exception &e) {
    std::cerr << "xml file " << xml_file << " contain errors" << std::endl;
  }
}

XmlDOMDocument::XmlDOMDocument(xercesc::DOMDocument* doc){
	doc_=doc;
}

XmlDOMDocument::~XmlDOMDocument() {
  if (doc_) {
    doc_->release();
  }
}

void XmlDOMDocument::addElement(xercesc::DOMElement* element){
	doc_->appendChild(element);
}

void XmlDOMDocument::addElement(xercesc::DOMElement* parent,xercesc::DOMElement* child){
	parent->appendChild(child);
}

int XmlDOMDocument::writeXmlDOMDocument(const char * filename){
	return prot::writeXmlFile(doc_,filename);
}

xercesc::DOMElement* getChildElement(xercesc::DOMElement *parent, 
                                const char* tag, int index) { 
  xercesc::DOMNodeList* list = parent->getElementsByTagName(X(tag));
  xercesc::DOMElement* element = 
      dynamic_cast<xercesc::DOMElement*>(list->item(index));
  if (element == nullptr) {
    std::stringstream stream;
    stream << "Get Child Element " << tag << " return null";
    LOG_ERROR(stream.str());
    throw stream.str();
  }
  return element;
}


std::string getChildValue(xercesc::DOMElement* parent,  
                          const char* child_tag, int i) {
  xercesc::DOMNodeList* node_list;
  node_list= parent->getElementsByTagName(X(child_tag));
  if (node_list == nullptr) {
    std::stringstream stream;
    stream << "Get Child Element " << child_tag << " return null";
    LOG_ERROR( stream.str());
    throw stream.str();
  }
  xercesc::DOMElement* child = 
  dynamic_cast<xercesc::DOMElement*>(node_list->item(i));
  if (child == nullptr) {
    std::stringstream stream;
    stream << "Get Child Element " << child_tag << " return null";
    LOG_ERROR( stream.str());
    throw stream.str();
  }

  std::string value;
  if (child) {
    value = Y(child->getTextContent());
  }
  else {
    value = "";
  }
  return value;
}

double getDoubleChildValue(xercesc::DOMElement* parent,  
                           const char* child_tag, int i) {
  std::string value = getChildValue(parent, child_tag, i);
  return atof(value.c_str());
}

int getIntChildValue(xercesc::DOMElement* parent,  
                     const char* child_tag, int i) {
  std::string value = getChildValue(parent, child_tag, i);
  return atoi(value.c_str());
}

bool getBoolChildValue(xercesc::DOMElement* parent,
                       const char* child_tag, int i) {
  std::string value = getChildValue(parent, child_tag, i);
  std::transform(value.begin(), value.end(), value.begin(), ::tolower);
  if(value == "true") {
    return true;
  }
  return false;
}

int getChildCount(xercesc::DOMElement* parent, 
                 const char* child_tag) {
  xercesc::DOMNodeList* childList = parent->getElementsByTagName(X(child_tag));
  return (int)childList->getLength();
}

std::string getAttributeValue(xercesc::DOMElement* element, const char* attribute_tag) {
  std::string value = Y(element->getAttribute(X(attribute_tag)));
  return value;
}


/**
 * Add an element 
 **/
xercesc::DOMElement* XmlDOMDocument::createElement(const char* tag) {
  xercesc::DOMElement* element = doc_->createElement(X(tag));
  return element;
}

xercesc::DOMText* XmlDOMDocument::createTextNode(const char* text) {
  xercesc::DOMText* text_node = doc_->createTextNode(X(text));
  return text_node;
}

void XmlDOMDocument::addElement(xercesc::DOMElement* element, 
                const char* tag, const char* value) {
  xercesc::DOMElement* child = createElement(tag);
  element->appendChild(child);
  xercesc::DOMText* text_node = createTextNode(value);
  child->appendChild(text_node);
}

std::string convertToString(double value) {
  std::stringstream stream;
  stream << value;
  return stream.str();
}

std::string convertToString(int value){
	std::stringstream stream;
	stream << value;
	return stream.str();
}

std::string convertToString(bool value) {
  std::stringstream stream;
  stream << value;
  return stream.str();
}

}

