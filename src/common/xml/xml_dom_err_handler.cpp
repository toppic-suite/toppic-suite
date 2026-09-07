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

#include "common/xml/xml_dom_err_handler.hpp"

#include <cstddef>
#include <pugixml.hpp>
#include <stdexcept>
#include <string>

#include "common/util/logger.hpp"

namespace toppic {

namespace xml_dom_err_handler {

namespace {

// Convert a byte offset in `source` to a 1-based "line L, column C" string.
// Falls back to the raw offset when `source` is empty or the offset is out of
// range (e.g. when the document was parsed straight from a file).
std::string describeLocation(const std::string& source, std::ptrdiff_t offset) {
  if (source.empty() || offset < 0 ||
      static_cast<std::size_t>(offset) > source.size()) {
    return "offset " + std::to_string(offset);
  }
  std::size_t line = 1;
  std::size_t column = 1;
  for (std::ptrdiff_t i = 0; i < offset; ++i) {
    if (source[static_cast<std::size_t>(i)] == '\n') {
      ++line;
      column = 1;
    } else {
      ++column;
    }
  }
  return "line " + std::to_string(line) + ", column " + std::to_string(column);
}

}  // namespace

void checkParseResult(const pugi::xml_parse_result& result,
                      const std::string& file_name, const std::string& source) {
  // xml_parse_result is contextually true only when status == status_ok.
  if (result) {
    return;
  }
  std::string location = describeLocation(source, result.offset);
  std::string message = "Fatal XML parsing error in " + file_name + " at " +
                        location + ": " + result.description();
  LOG_ERROR(message);
  throw std::runtime_error(message);
}

}  // namespace xml_dom_err_handler

}  // namespace toppic
