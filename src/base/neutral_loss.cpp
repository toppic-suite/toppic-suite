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

#include "xml/xml_dom_util.hpp"
#include "base/neutral_loss.hpp"

namespace toppic {

NeutralLoss::NeutralLoss(const std::string &name, double mass): 
    name_(name), mass_(mass) { }

NeutralLoss::NeutralLoss(XmlDOMElement* element) {
  name_ = xml_dom_util::getChildValue(element, "name", 0);
  mass_ = xml_dom_util::getDoubleChildValue(element, "mass", 0);
}

}  // namespace toppic
