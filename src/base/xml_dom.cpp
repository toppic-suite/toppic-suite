//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#include <iostream>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/util/OutOfMemoryException.hpp>
#include <xercesc/framework/XMLFormatter.hpp>
#include <xercesc/framework/LocalFileFormatTarget.hpp>
#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMImplementation.hpp>
#include <xercesc/dom/DOMImplementationRegistry.hpp>
#include <xercesc/dom/DOMLSSerializer.hpp>
#include <xercesc/dom/DOMLSOutput.hpp>

#include "base/xml_dom.hpp"
#include "base/xml_dom_err_handler.hpp"
 
namespace prot {

XmlDOMParser* XmlDOMParserFactory::dom_parser_ = nullptr;


/* XmlDOMParser */
XmlDOMParser::XmlDOMParser() : parser_(nullptr), err_handler_(nullptr) {
  xercesc::XMLPlatformUtils::Initialize();
  parser_ = new xercesc::XercesDOMParser();
  err_handler_ = (xercesc::ErrorHandler*) new XmlDOMErrorHandler();
  parser_->setErrorHandler(err_handler_);
}

XmlDOMParser::~XmlDOMParser() {
  if (parser_ != nullptr) {
    delete parser_;
    xercesc::XMLPlatformUtils::Terminate();
  }
}

xercesc::DOMDocument* XmlDOMParser::parse(const std::string &xml_file) {
  parser_->parse(xml_file.c_str());
  return parser_->adoptDocument();
}

xercesc::DOMDocument* XmlDOMParser::parse(const xercesc::MemBufInputSource &str_buf) {
  parser_->parse(str_buf);
  return parser_->adoptDocument();
}

/* XmlDOMParserFactory */
XmlDOMParser* XmlDOMParserFactory::getXmlDOMParserInstance() {
  if (dom_parser_ == nullptr) {
    dom_parser_ = new XmlDOMParser();
  }
  return dom_parser_;
}

/* XmlDOMImplenmation */
XmlDOMImpl* XmlDOMImplFactory::dom_impl_ = nullptr;

XmlDOMImpl::XmlDOMImpl() {
  impl_ = xercesc::DOMImplementationRegistry::getDOMImplementation(X("Core"));
}

XmlDOMImpl::~XmlDOMImpl() {
  if (impl_ != nullptr) {
    delete impl_;
  }
}

xercesc::DOMDocument* XmlDOMImpl::createDoc(const std::string &root) {
  xercesc::DOMDocument* doc = impl_->createDocument(0, X(root.c_str()), 0);
  return doc;
}

xercesc::DOMLSSerializer* XmlDOMImpl::createSerializer() {
  xercesc::DOMLSSerializer* writer = impl_->createLSSerializer();
  writer->getDomConfig()->setParameter(
      xercesc::XMLUni::fgDOMWRTFormatPrettyPrint, true);
  writer->getDomConfig()->setParameter(
      xercesc::XMLUni::fgDOMWRTDiscardDefaultContent, true);
  writer->setNewLine(X("\n"));
  return writer;
}

/*XmlDOMImplFactory */
XmlDOMImpl* XmlDOMImplFactory::getXmlDOMImplInstance() {
  if (dom_impl_ == nullptr) {
    dom_impl_ = new XmlDOMImpl();
  }
  return dom_impl_;
}

}
