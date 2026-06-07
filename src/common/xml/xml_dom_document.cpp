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

#include "common/xml/xml_dom_document.hpp"

#include <memory>
#include <string>
#include <utility>

#include <pugixml.hpp>

namespace toppic {

XmlDOMDocument::XmlDOMDocument(std::unique_ptr<pugi::xml_document> doc)
    : doc_(std::move(doc)) {}

XmlDOMDocument::XmlDOMDocument(const std::string &root)
    : doc_(std::make_unique<pugi::xml_document>()) {
  doc_->append_child(root.c_str());
}

XmlDOMElement XmlDOMDocument::getDocumentElement() {
  return doc_->document_element();
}

XmlDOMElement XmlDOMDocument::addElement(XmlDOMElement parent, const char* tag) {
  return parent.append_child(tag);
}

void XmlDOMDocument::addElement(XmlDOMElement parent,
                                const char* tag, const char* value) {
  parent.append_child(tag).text().set(value);
}

}  // namespace toppic
