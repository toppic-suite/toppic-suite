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

#ifndef TOPPIC_COMMON_XML_XML_DOM_DOCUMENT_HPP_
#define TOPPIC_COMMON_XML_XML_DOM_DOCUMENT_HPP_

#include <memory>
#include <string>

#include <pugixml.hpp>

#include "common/xml/xml_dom_element.hpp"

namespace toppic {

// Owns a pugi::xml_document and exposes element creation/lookup. Unlike the
// Xerces version this is RAII (no manual release()), and because pugixml has no
// detached nodes, createElement()/createTextNode() are replaced by an
// addElement() that creates the child already attached to its parent.
class XmlDOMDocument {
 public:
  // Take ownership of a parsed document (from XmlDOMParser::parse / parseStr).
  explicit XmlDOMDocument(std::unique_ptr<pugi::xml_document> doc);

  // Create a new document with a single root element named `root`.
  explicit XmlDOMDocument(const std::string &root);

  XmlDOMDocument(const XmlDOMDocument&) = delete;
  XmlDOMDocument& operator=(const XmlDOMDocument&) = delete;

  // The root element of the document.
  XmlDOMElement getDocumentElement();

  // Append a child element named `tag` to `parent` and return it (build its
  // subtree by passing the returned node as the next parent). pugixml creates
  // the node already attached, so there is no separate "append" step.
  XmlDOMElement addElement(XmlDOMElement parent, const char* tag);

  // Append a leaf element <tag>value</tag> to `parent`.
  void addElement(XmlDOMElement parent, const char* tag, const char* value);

 private:
  std::unique_ptr<pugi::xml_document> doc_;
};

}  // namespace toppic

#endif  // TOPPIC_COMMON_XML_XML_DOM_DOCUMENT_HPP_
