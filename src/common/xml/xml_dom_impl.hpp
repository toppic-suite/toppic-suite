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

#ifndef TOPPIC_COMMON_XML_XML_DOM_IMPL_HPP_
#define TOPPIC_COMMON_XML_XML_DOM_IMPL_HPP_

#include <memory>
#include <string>

#include <pugixml.hpp>

namespace toppic {

// Xerces needed a DOMImplementation singleton (the old XmlDOMImpl /
// XmlDOMImplFactory) to create documents and serializers. pugixml needs
// neither: a pugi::xml_document is a self-contained value, and serialization is
// a node method (see xml_dom_util::writeToString). So the serializer factory
// and the singleton are gone; the only piece still worth a helper is creating a
// document with a named root.
namespace xml_dom_impl {

// Returns a new document whose root element is named `root`.
std::unique_ptr<pugi::xml_document> createDoc(const std::string &root);

}  // namespace xml_dom_impl

}  // namespace toppic

#endif  // TOPPIC_COMMON_XML_XML_DOM_IMPL_HPP_
