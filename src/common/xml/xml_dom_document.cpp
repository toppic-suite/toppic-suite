//Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane University.
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

#include "common/xml/xml_dom_document.hpp"

#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMText.hpp>

#include "common/util/logger.hpp"
#include "common/xml/xml_dom_str.hpp"

namespace toppic {

XmlDOMDocument::XmlDOMDocument(XmlDOMParser* parser,
                               const char* xml_file)
    : doc_(nullptr) {
  try {
    doc_ = parser->parse(xml_file);
  }
  catch (std::exception &e) {
    LOG_ERROR("The xml file " << xml_file << " contains errors: " << e.what());
    throw;
  }
}

XmlDOMDocument::XmlDOMDocument(XmlDOMParser* parser,
                               const xercesc::MemBufInputSource &str_buf)
    : doc_(nullptr) {
  try {
    doc_ = parser->parse(str_buf);
  }
  catch (std::exception &e) {
    LOG_ERROR("Xml str buffer contains errors: " << e.what());
    throw;
  }
}

XmlDOMDocument::XmlDOMDocument(xercesc::DOMImplementation* implementation,
                               const std::string &root) {
  doc_ = implementation->createDocument(0, XmlStr(root.c_str()).unicodeForm(), 0);
}

XmlDOMDocument::XmlDOMDocument(xercesc::DOMDocument* doc) {
  doc_ = doc;
}

XmlDOMDocument::~XmlDOMDocument() {
  if (doc_) {
    doc_->release();
  }
}

void XmlDOMDocument::addElement(XmlDOMElement* element) {
  doc_->appendChild(element);
}

void XmlDOMDocument::addElement(XmlDOMElement* parent, XmlDOMElement* child) {
  parent->appendChild(child);
}

XmlDOMElement* XmlDOMDocument::createElement(const char* tag) {
  return doc_->createElement(XmlStr(tag).unicodeForm());
}

xercesc::DOMText* XmlDOMDocument::createTextNode(const char* text) {
  return doc_->createTextNode(XmlStr(text).unicodeForm());
}

void XmlDOMDocument::addElement(XmlDOMElement* element,
                                const char* tag, const char* value) {
  XmlDOMElement* child = createElement(tag);
  element->appendChild(child);
  xercesc::DOMText* text_node = createTextNode(value);
  child->appendChild(text_node);
}

}  // namespace toppic
