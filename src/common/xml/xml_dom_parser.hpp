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

#ifndef TOPPIC_COMMON_XML_XML_DOM_PARSER_HPP_
#define TOPPIC_COMMON_XML_XML_DOM_PARSER_HPP_

#include <memory>
#include <string>

#include <pugixml.hpp>

namespace toppic {

// pugixml documents parse themselves (load_file / load_buffer) and need no
// global initialization (Xerces required XMLPlatformUtils + an attached
// ErrorHandler). So this is a thin wrapper that adds the project's parse-error
// handling (xml_dom_err_handler). The factory/singleton is kept only for
// source-compatibility with the existing call sites.
class XmlDOMParser {
 public:
  XmlDOMParser() = default;
  XmlDOMParser(const XmlDOMParser&) = delete;
  XmlDOMParser& operator=(const XmlDOMParser&) = delete;

  // Parse an XML file; throws std::runtime_error on a parse error.
  std::unique_ptr<pugi::xml_document> parse(const std::string &xml_file);

  // Parse an in-memory XML string; throws std::runtime_error on a parse error.
  std::unique_ptr<pugi::xml_document> parseStr(const std::string &xml_str);
};

class XmlDOMParserFactory {
 public:
  XmlDOMParserFactory() = delete;
  static XmlDOMParser* getXmlDOMParserInstance();
  static void deleteParserInstance();

 private:
  static XmlDOMParser* dom_parser_;
};

}  // namespace toppic

#endif  // TOPPIC_COMMON_XML_XML_DOM_PARSER_HPP_
