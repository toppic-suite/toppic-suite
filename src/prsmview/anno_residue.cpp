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


#include "base/string_util.hpp"
#include "prsmview/anno_residue.hpp"

namespace prot {

AnnoResidue::AnnoResidue(ResiduePtr residue_ptr, int pos):
    Residue(residue_ptr->getAcidPtr(),residue_ptr->getPtmPtr()) {
      pos_= pos;
      type_= ANNO_RESIDUE_TYPE_NORMAL;
      is_unexpected_change_ = false;
      unexpected_change_color_ = 0;
    }

void AnnoResidue::appendViewXml(XmlDOMDocument* xml_doc,
                                xercesc::DOMElement* parent){
  xercesc::DOMElement* element = xml_doc->createElement("residue");
  std::string str = string_util::convertToString(pos_);
  xml_doc->addElement(element, "position", str.c_str());
  str = getAcidPtr()->getOneLetter();
  xml_doc->addElement(element, "acid", str.c_str());
  str = type_;
  xml_doc->addElement(element, "residue_type", str.c_str());
  str = string_util::convertToString(is_unexpected_change_);
  xml_doc->addElement(element, "is_unexpected_change", str.c_str());
  str = string_util::convertToString(unexpected_change_color_);
  xml_doc->addElement(element, "unexpected_change_color", str.c_str());
  str = string_util::convertToString(possible_pos_color_);
  xml_doc->addElement(element, "possible_pos_color", str.c_str());
  xml_doc->addElement(element, "anno", anno_.c_str());
  parent->appendChild(element);
}

}

