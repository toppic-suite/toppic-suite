//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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
#include "prsmview/anno_ptm_position.hpp"

namespace toppic {

AnnoPtmPosition::AnnoPtmPosition(int left_pos, int right_pos, 
                                 std::string anno): 
    left_pos_(left_pos), 
    right_pos_(right_pos), 
    anno_(anno) {}

void AnnoPtmPosition::appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent) {
  xercesc::DOMElement* position_element = xml_doc->createElement("occurence");
  std::string str = str_util::toString(left_pos_);
  xml_doc->addElement(position_element, "left_pos", str.c_str());
  str = str_util::toString(right_pos_);
  xml_doc->addElement(position_element, "right_pos", str.c_str());
  xml_doc->addElement(position_element, "anno", anno_.c_str());

  parent->appendChild(position_element);
}

}  // namespace toppic
