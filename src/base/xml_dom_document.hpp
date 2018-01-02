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


#ifndef PROT_BASE_XML_DOM_DOCUMENT_HPP_
#define PROT_BASE_XML_DOM_DOCUMENT_HPP_

#include <iostream>
#include <fstream>
#include <string>

#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/sax/HandlerBase.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/framework/MemBufInputSource.hpp>

#include "base/xml_dom.hpp"

namespace prot {

class XmlDOMDocument {

 public:
  XmlDOMDocument(XmlDOMParser* parser, const char* xml_file);
  XmlDOMDocument(XmlDOMParser* parser, const xercesc::MemBufInputSource &str_buf);

  XmlDOMDocument(xercesc::DOMDocument* doc);

  XmlDOMDocument(xercesc::DOMImplementation* implementation, 
                 const std::string &root);
  ~XmlDOMDocument();

  xercesc::DOMElement* createElement(const char* tag);

  xercesc::DOMText* createTextNode(const char* text);

  void addElement(xercesc::DOMElement* element, 
                  const char* tag, const char* value);

  xercesc::DOMElement* getDocumentElement() {
    return doc_->getDocumentElement();
  }

  void addElement(xercesc::DOMElement* element);
  void addElement(xercesc::DOMElement* parent,xercesc::DOMElement* child);

  int writeXmlDOMDocument(const char * filename);

 private:
  xercesc::DOMDocument* doc_;
};

}

#endif
