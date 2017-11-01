//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#include "base/mass_constant.hpp"
#include "base/ion_type.hpp"
#include "base/string_util.hpp"
#include "base/xml_dom_document.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

IonType::IonType(const std::string &name, bool n_term, double shift): 
    name_(name),
    n_term_(n_term),
    shift_(shift) {
      if (n_term_) {
        b_y_shift_ = shift_;
      }
      else {
        b_y_shift_ = shift_ - MassConstant::getYIonShift();
      }
    }

IonType::IonType(xercesc::DOMElement* element) { 
  name_ = XmlDomUtil::getChildValue(element, "name", 0);
  n_term_ = XmlDomUtil::getBoolChildValue(element, "n_term", 0);
  shift_ = XmlDomUtil::getDoubleChildValue(element, "shift", 0);
  if (n_term_) {
    b_y_shift_ = shift_;
  }
  else {
    b_y_shift_ = shift_ - MassConstant::getYIonShift();
  }
}

void IonType::appendNameToXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
  xercesc::DOMElement* element = xml_doc->createElement("ion_type");
  xml_doc->addElement(element, "name", name_.c_str());
  parent->appendChild(element);
}

}
