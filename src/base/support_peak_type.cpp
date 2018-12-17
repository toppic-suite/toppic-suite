//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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


#include "base/support_peak_type.hpp"
#include "base/xml_dom_util.hpp"

namespace toppic {

SupportPeakType::SupportPeakType(xercesc::DOMElement* element) {
  id_ = xml_dom_util::getIntChildValue(element, "id", 0);
  name_ = xml_dom_util::getChildValue(element, "name", 0);
}

}  // namespace toppic
