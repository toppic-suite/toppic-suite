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

}
