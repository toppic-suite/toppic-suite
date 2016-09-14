// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


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
