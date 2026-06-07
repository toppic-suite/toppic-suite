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

#ifndef TOPPIC_COMMON_XML_XML_DOM_ERR_HANDLER_HPP_
#define TOPPIC_COMMON_XML_XML_DOM_ERR_HANDLER_HPP_

#include <string>

#include <pugixml.hpp>

namespace toppic {

namespace xml_dom_err_handler {

// pugixml reports parse problems through xml_parse_result (a status code, a
// human-readable description and a byte offset) rather than a SAX error-handler
// callback, so this replaces the old Xerces XmlDOMErrorHandler::fatalError.
//
// On a failed parse it logs the error and throws std::runtime_error. file_name
// is used only for the message. If `source` (the parsed text) is given, the
// byte offset is reported as a 1-based line:column; otherwise the raw offset is
// reported.
void checkParseResult(const pugi::xml_parse_result& result,
                      const std::string& file_name,
                      const std::string& source = std::string());

}  // namespace xml_dom_err_handler

}  // namespace toppic

#endif  // TOPPIC_COMMON_XML_XML_DOM_ERR_HANDLER_HPP_
