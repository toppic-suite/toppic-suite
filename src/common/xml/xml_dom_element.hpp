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

#ifndef TOPPIC_COMMON_XML_XML_DOM_ELEMENT_HPP_
#define TOPPIC_COMMON_XML_XML_DOM_ELEMENT_HPP_

#include <pugixml.hpp>

namespace toppic {

// pugixml uses lightweight value handles (pugi::xml_node) rather than DOM node
// pointers, so XmlDOMElement is now a value type. Code that previously passed
// XmlDOMElement* now passes XmlDOMElement (by value or const reference); an
// empty/absent node is tested with `if (!element)` instead of `== nullptr`.
using XmlDOMElement = pugi::xml_node;

}  // namespace toppic

#endif
