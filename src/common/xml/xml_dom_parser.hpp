//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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

#ifndef TOPPIC_COMMON_XML_XML_DOM_PARSER_HPP_
#define TOPPIC_COMMON_XML_XML_DOM_PARSER_HPP_

#include <string>

#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/sax/HandlerBase.hpp>
#include <xercesc/framework/MemBufInputSource.hpp>

namespace toppic {

// DOM parser
class XmlDOMParser {
 public:
  XmlDOMParser();
  ~XmlDOMParser();

  xercesc::DOMDocument* parse(const std::string &xml_file);

  xercesc::DOMDocument* parse(const xercesc::MemBufInputSource &str_buf);

 private:
  xercesc::XercesDOMParser* parser_;
  xercesc::ErrorHandler*    err_handler_;
};

class XmlDOMParserFactory {
 private:
  static XmlDOMParser* dom_parser_;
 public:
  static XmlDOMParser* getXmlDOMParserInstance();
};

}  // namespace toppic

#endif
