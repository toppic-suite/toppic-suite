#include <string>
#include <sstream>
#include <exception>
#include <algorithm>
#include <iomanip>


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

XmlDOMDocument::XmlDOMDocument(XmlDOMParser* parser, const xercesc::MemBufInputSource &str_buf) {
  try {
    doc_ = parser->parse(str_buf);
  }
  catch (std::exception &e) {
    std::cerr << "xml str buffer contain errors" << std::endl;
  }
}

XmlDOMDocument::XmlDOMDocument(xercesc::DOMImplementation* implementation,
                               const std::string &root){
    xercesc::DOMDocument* doc = implementation->createDocument(0,X(root.c_str()),0);
    doc_=doc;
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

xercesc::DOMNodeList* getChildElements(xercesc::DOMElement *parent,
                                const char* tag) {
  return parent->getElementsByTagName(X(tag));
}

xercesc::DOMElement* getChildElement(xercesc::DOMElement *parent, 
                                const char* tag, int index) { 
  xercesc::DOMNodeList* list = parent->getElementsByTagName(X(tag));
  xercesc::DOMElement* element = 
      dynamic_cast<xercesc::DOMElement*>(list->item(index));
  if (element == nullptr) {
    std::stringstream stream;
    stream << "Get Child Element " << tag << " return null";
    LOG_WARN(stream.str());
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
    LOG_WARN( stream.str());
    throw stream.str();
  }
  xercesc::DOMElement* child = 
  dynamic_cast<xercesc::DOMElement*>(node_list->item(i));
  if (child == nullptr) {
    std::stringstream stream;
    stream << "Get Child Element " << child_tag << " return null";
    LOG_WARN( stream.str());
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
  return std::stod(value);
}

int getIntChildValue(xercesc::DOMElement* parent,  
                     const char* child_tag, int i) {
  try {
    std::string value = getChildValue(parent, child_tag, i);
    return std::stoi(value);
  }
  catch (std::string s) {
    return 0;
  }
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

  if(value < 1 && value > -1 && value !=0){
    stream << std::scientific << std::setprecision(10);
  }
  else{
    stream << std::fixed<<std::setprecision(10);
  }
  stream << value;
  return stream.str();
}

std::string convertToString(double value, int number) {
  std::stringstream stream;
  if(value ==0)
  {
    stream << std::fixed << std::setprecision(0);
  }
  else if (value < 0.01 && value > -0.01 && value != 0) {
    if(number>2){
      stream << std::scientific << std::setprecision(2);
    }
    else{
      stream << std::scientific << std::setprecision(number);
    }
  } else {
    stream << std::fixed << std::setprecision(number);
  }
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

std::string writeToString(xercesc::DOMLSSerializer* serializer, xercesc::DOMNode *node) {
  XMLCh* ch = serializer->writeToString(node, 0);
  std::string result = Y(ch);
  xercesc::XMLString::release(&ch);
  return result;
}

void writeToStreamByRemovingDoubleLF(std::ofstream &file, std::string &str) {
  int pos = 0;
  std::size_t found = str.find("\n\n", pos);
  while (found != std::string::npos) {
    std::string sub = str.substr(pos, found - pos);
    file << sub << std::endl;
    pos = found + 2;
    found = str.find("\n\n", pos);
  }
  if (pos < (int)str.length()) {
    std::string sub = str.substr(pos);
    file << sub << std::endl; 
  }
}

}
