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

#include "common/xml/xml_dom_impl.hpp"

#include <memory>
#include <string>

#include <pugixml.hpp>

namespace toppic {

namespace xml_dom_impl {

std::unique_ptr<pugi::xml_document> createDoc(const std::string &root) {
  auto doc = std::make_unique<pugi::xml_document>();
  doc->append_child(root.c_str());
  return doc;
}

}  // namespace xml_dom_impl

}  // namespace toppic
