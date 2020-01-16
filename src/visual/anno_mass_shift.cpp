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

#include "common/util/str_util.hpp"
#include "visual/anno_mass_shift.hpp"

namespace toppic {

AnnoMassShift::AnnoMassShift(int id, int left_pos, int right_pos, 
                             const std::string & anno_str, 
                             AlterTypePtr & mass_shift_type):
    id_(id), 
    left_pos_(left_pos),
    right_pos_(right_pos),
    anno_str_(anno_str),
    mass_shift_type_(mass_shift_type) {}

void AnnoMassShift::appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent) {
  xercesc::DOMElement* element = xml_doc->createElement("mass_shift");
  std::string str = str_util::toString(id_);
  xml_doc->addElement(element, "id", str.c_str());
  str = str_util::toString(left_pos_);
  xml_doc->addElement(element, "left_position", str.c_str());
  str = str_util::toString(right_pos_);
  xml_doc->addElement(element, "right_position", str.c_str());

  xml_doc->addElement(element, "anno", anno_str_.c_str());

  std::string type; 
  if (mass_shift_type_ == toppic::AlterType::UNEXPECTED) {
    type = "unexpected";
  }
  else {
    type = "variable ptm";
  }
  xml_doc->addElement(element, "shift_type", type.c_str());
  parent->appendChild(element);
}

}  // namespace toppic

