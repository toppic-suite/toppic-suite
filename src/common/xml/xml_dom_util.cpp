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

#include "common/xml/xml_dom_util.hpp"

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <sstream>
#include <stdexcept>
#include <string>

#include <pugixml.hpp>

#include "common/util/logger.hpp"

namespace toppic {

namespace xml_dom_util {

XmlDOMElement getChildElement(const XmlDOMElement& parent,
                              const char* tag, int index) {
  // Direct children named `tag`, in document order. The call-site audit
  // confirmed no query relies on a deeper (descendant) match, so this is a
  // faster, behavior-equivalent replacement for an XPath ".//tag" search.
  int count = 0;
  for (pugi::xml_node child : parent.children(tag)) {
    if (count == index) {
      return child;
    }
    ++count;
  }
  LOG_WARN("Get Child Element " << tag << " return null!");
  throw std::runtime_error(std::string("getChildElement: element not found: ") + tag);
}

std::string getChildValue(const XmlDOMElement& parent,
                          const char* child_tag, int i) {
  XmlDOMElement child = getChildElement(parent, child_tag, i);
  // text() returns the element's character data (e.g. "foo" in <t>foo</t>).
  return child.text().as_string();
}

double getScientificChildValue(const XmlDOMElement& parent,
                               const char* child_tag, int i) {
  std::string value = getChildValue(parent, child_tag, i);
  return std::stod(value);
}

double getDoubleChildValue(const XmlDOMElement& parent,
                           const char* child_tag, int i) {
  std::string value = getChildValue(parent, child_tag, i);
  return std::stod(value);
}

int getIntChildValue(const XmlDOMElement& parent,
                     const char* child_tag, int i) {
  try {
    std::string value = getChildValue(parent, child_tag, i);
    return std::stoi(value);
  }
  catch (const std::logic_error& e) {
    LOG_WARN("Get Child Element " << child_tag << " error: " << e.what());
    return 0;
  }
}

bool getBoolChildValue(const XmlDOMElement& parent,
                       const char* child_tag, int i) {
  std::string value = getChildValue(parent, child_tag, i);
  // Cast to unsigned char: passing a negative char to std::tolower is UB.
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return value == "true";
}

int getChildCount(const XmlDOMElement& parent, const char* child_tag) {
  int count = 0;
  for (pugi::xml_node child : parent.children(child_tag)) {
    (void) child;
    ++count;
  }
  return count;
}

std::string getAttributeValue(const XmlDOMElement& element,
                              const char* attribute_tag) {
  // as_string() returns "" when the attribute is absent (as Xerces did).
  return element.attribute(attribute_tag).as_string();
}

std::string writeToString(const XmlDOMElement& node) {
  std::ostringstream stream;
  node.print(stream, "  ");
  return stream.str();
}

void writeToStreamByRemovingDoubleLF(std::ofstream& file, const std::string& str) {
  std::size_t pos = 0;
  std::size_t found = str.find("\n\n", pos);
  while (found != std::string::npos) {
    std::string sub = str.substr(pos, found - pos);
    file << sub << std::endl;
    pos = found + 2;
    found = str.find("\n\n", pos);
  }
  if (pos < str.length()) {
    std::string sub = str.substr(pos);
    file << sub << std::endl;
  }
}

}  // namespace xml_dom_util

}  // namespace toppic
