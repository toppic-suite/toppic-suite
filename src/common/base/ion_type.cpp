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

#include <string>

#include "common/util/str_util.hpp"
#include "common/xml/xml_dom_document.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "common/base/mass_constant.hpp"
#include "common/base/ion_type.hpp"

namespace toppic {

IonType::IonType(const std::string &name, bool n_term, double shift):
    name_(name),
    n_term_(n_term),
    shift_(shift) {
      if (n_term_) {
        b_y_shift_ = shift_;
      } else {
        b_y_shift_ = shift_ - mass_constant::getYIonShift();
      }
    }

IonType::IonType(XmlDOMElement* element) {
  name_ = xml_dom_util::getChildValue(element, "name", 0);
  n_term_ = xml_dom_util::getBoolChildValue(element, "n_term", 0);
  shift_ = xml_dom_util::getDoubleChildValue(element, "shift", 0);
  if (n_term_) {
    b_y_shift_ = shift_;
  } else {
    b_y_shift_ = shift_ - mass_constant::getYIonShift();
  }
}

void IonType::appendNameToXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent) {
  XmlDOMElement* element = xml_doc->createElement("ion_type");
  xml_doc->addElement(element, "name", name_.c_str());
  parent->appendChild(element);
}

}  // namespace toppic
