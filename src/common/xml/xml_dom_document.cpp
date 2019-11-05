//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMText.hpp>

#include "common/util/logger.hpp"
#include "common/xml/xml_dom_str.hpp"
#include "common/xml/xml_dom_document.hpp"

namespace toppic {

XmlDOMDocument::XmlDOMDocument(XmlDOMParser* parser, 
                               const char* xml_file) : doc_(NULL) {
  try {
    doc_ = parser->parse(xml_file);
  }
  catch (std::exception &e) {
    LOG_ERROR("xml file " << xml_file << " contain errors!");
    exit(EXIT_FAILURE);
  }
}

XmlDOMDocument::XmlDOMDocument(XmlDOMParser* parser, const xercesc::MemBufInputSource &str_buf) {
  try {
    doc_ = parser->parse(str_buf);
  }
  catch (std::exception &e) {
    LOG_ERROR("xml str buffer contain errors!");
    exit(EXIT_FAILURE);
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
