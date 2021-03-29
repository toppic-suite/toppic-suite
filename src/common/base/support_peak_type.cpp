//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
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
#include "common/base/support_peak_type.hpp"

namespace toppic {

SupportPeakType::SupportPeakType(int id, const std::string &name): 
    id_(id), 
    name_(name) {}

SupportPeakType::SupportPeakType(XmlDOMElement* element) {
  id_ = xml_dom_util::getIntChildValue(element, "id", 0);
  name_ = xml_dom_util::getChildValue(element, "name", 0);
}

}  // namespace toppic
