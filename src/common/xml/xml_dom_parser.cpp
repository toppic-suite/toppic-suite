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

#include "common/xml/xml_dom_parser.hpp"

#include <xercesc/util/PlatformUtils.hpp>

#include "common/xml/xml_dom_err_handler.hpp"

namespace toppic {

XmlDOMParser* XmlDOMParserFactory::dom_parser_ = nullptr;

/* XmlDOMParser */
XmlDOMParser::XmlDOMParser() : parser_(nullptr), err_handler_(nullptr) {
  xercesc::XMLPlatformUtils::Initialize();
  try {
    parser_ = new xercesc::XercesDOMParser();
    err_handler_ = new XmlDOMErrorHandler();
    parser_->setErrorHandler(err_handler_);
  } catch (...) {
    delete err_handler_;
    delete parser_;
    xercesc::XMLPlatformUtils::Terminate();
    throw;
  }
}

XmlDOMParser::~XmlDOMParser() {
  delete parser_;
  delete err_handler_;
  xercesc::XMLPlatformUtils::Terminate();
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

void XmlDOMParserFactory::deleteParserInstance() {
  delete dom_parser_;
  dom_parser_ = nullptr;
}

}  // namespace toppic
