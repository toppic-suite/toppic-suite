// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "common/xml/xml_dom_parser.hpp"

#include <memory>
#include <pugixml.hpp>
#include <string>

#include "common/xml/xml_dom_err_handler.hpp"

namespace toppic {

XmlDOMParser* XmlDOMParserFactory::dom_parser_ = nullptr;

std::unique_ptr<pugi::xml_document> XmlDOMParser::parse(
    const std::string& xml_file) {
  auto doc = std::make_unique<pugi::xml_document>();
  pugi::xml_parse_result result = doc->load_file(xml_file.c_str());
  // The file is read by pugixml, so the source text is not available here; the
  // error location is reported as a byte offset.
  xml_dom_err_handler::checkParseResult(result, xml_file);
  return doc;
}

std::unique_ptr<pugi::xml_document> XmlDOMParser::parseStr(
    const std::string& xml_str) {
  auto doc = std::make_unique<pugi::xml_document>();
  pugi::xml_parse_result result =
      doc->load_buffer(xml_str.data(), xml_str.size());
  xml_dom_err_handler::checkParseResult(result, "<in-memory XML>", xml_str);
  return doc;
}

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
