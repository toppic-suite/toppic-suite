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
#include <memory>
#include <stdexcept>
#include <string>

#include <xercesc/dom/DOMLSSerializer.hpp>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include <xercesc/util/XMLString.hpp>

#include "common/util/logger.hpp"
#include "common/xml/xml_dom_str.hpp"

namespace toppic {

namespace xml_dom_util {

XmlDOMElement* getChildElement(XmlDOMElement* parent,
                               const char* tag, int index) {
  xercesc::DOMNodeList* list = parent->getElementsByTagName(XmlStr(tag).unicodeForm());
  XmlDOMElement* element = dynamic_cast<XmlDOMElement*>(list->item(index));
  if (element == nullptr) {
    LOG_WARN("Get Child Element " << tag << " return null!");
    throw std::runtime_error(std::string("getChildElement: element not found: ") + tag);
  }
  return element;
}

std::string getChildValue(XmlDOMElement* parent,
                          const char* child_tag, int i) {
  xercesc::DOMNodeList* node_list = parent->getElementsByTagName(XmlStr(child_tag).unicodeForm());
  if (node_list->getLength() == 0) {
    LOG_WARN("Get Child Element " << child_tag << " return null!");
    throw std::runtime_error(std::string("getChildValue: node list not found: ") + child_tag);
  }
  XmlDOMElement* child = dynamic_cast<XmlDOMElement*>(node_list->item(i));
  if (child == nullptr) {
    LOG_WARN("Get Child Element " << child_tag << " return null!");
    throw std::runtime_error(std::string("getChildValue: element not found: ") + child_tag);
  }
  return CharStr(child->getTextContent()).getString();
}

double getScientificChildValue(XmlDOMElement* parent,
                               const char* child_tag, int i) {
  std::string value = getChildValue(parent, child_tag, i);
  // std::stod parses scientific notation (e.g. "1.5e-3"), so it is a drop-in
  // for the former str_util::scientificToDouble.
  return std::stod(value);
}

double getDoubleChildValue(XmlDOMElement* parent,
                           const char* child_tag, int i) {
  std::string value = getChildValue(parent, child_tag, i);
  return std::stod(value);
}

int getIntChildValue(XmlDOMElement* parent,
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

bool getBoolChildValue(XmlDOMElement* parent,
                       const char* child_tag, int i) {
  std::string value = getChildValue(parent, child_tag, i);
  // Cast to unsigned char: passing a negative char to std::tolower is UB.
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return value == "true";
}

int getChildCount(XmlDOMElement* parent, const char* child_tag) {
  xercesc::DOMNodeList* child_list = parent->getElementsByTagName(XmlStr(child_tag).unicodeForm());
  return static_cast<int>(child_list->getLength());
}

std::string getAttributeValue(XmlDOMElement* element,
                              const char* attribute_tag) {
  return CharStr(element->getAttribute(XmlStr(attribute_tag).unicodeForm())).getString();
}

namespace {
// Deleter so an XMLCh* buffer allocated by Xerces (e.g. by
// DOMLSSerializer::writeToString) can be owned by a std::unique_ptr.
struct XmlChDeleter {
  void operator()(XMLCh* p) const { xercesc::XMLString::release(&p); }
};
}  // namespace

std::string writeToString(xercesc::DOMLSSerializer* serializer,
                          xercesc::DOMNode* node) {
  XMLCh* raw = serializer->writeToString(node, 0);
  if (raw == nullptr) {
    throw std::runtime_error("writeToString: serializer returned null");
  }
  // RAII: the buffer is released even if the transcoding below throws.
  std::unique_ptr<XMLCh, XmlChDeleter> ch(raw);
  return CharStr(ch.get()).getString();
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
