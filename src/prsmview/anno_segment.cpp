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

#include <cmath>

#include "base/string_util.hpp"
#include "prsmview/anno_segment.hpp"

namespace prot {

AnnoSegment::AnnoSegment(std::string segment_type, int left_pos, int right_pos,
                         double mass_shift, int color) {
  segment_type_ = segment_type;
  left_pos_ = left_pos;
  right_pos_ = right_pos;
  mass_shift_ = mass_shift;
  color_ = color;
}

void AnnoSegment::addOccurence(int pos, const std::string &acid_letter) {
  std::pair<int, std::string> new_occurence(pos + 1, acid_letter);
  occurences_.push_back(new_occurence);
}

std::string AnnoSegment::getResidueAnno() {
  if (anno_ != "") {
    return anno_;
  }

  occu_ = "";
  if (ptm_ptr_ != nullptr) {
    anno_ += "PTM: " + ptm_ptr_->getName() + "\n";
    for (size_t i = 0; i < occurences_.size(); i++) {
      if (score_[i] == 1) score_[i] = 0.999;

      if (score_[i] == 0) continue;

      score_[i] = std::floor(score_[i] * 1000) / 10;

      anno_ += "Site: " + occurences_[i].second + StringUtil::convertToString(occurences_[i].first) + " ";
      anno_ += "Confidence: " + StringUtil::convertToString(score_[i], 1) + "%\n";
      occu_ += occurences_[i].second + StringUtil::convertToString(occurences_[i].first) 
          + ":" + StringUtil::convertToString(score_[i], 1) + "%";
      if (i != occurences_.size() - 1) {
        occu_ += "; ";
      } 
    }
  }
  return anno_;
}

void AnnoSegment::appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent,
                            int precison){
  xercesc::DOMElement* element = xml_doc->createElement("unexpected_change");
  std::string str = StringUtil::convertToString(left_pos_);
  xml_doc->addElement(element, "left_position", str.c_str());
  str = StringUtil::convertToString(right_pos_);
  xml_doc->addElement(element, "right_position", str.c_str());
  str = StringUtil::convertToString(mass_shift_, precison);
  xml_doc->addElement(element, "mass_shift", str.c_str());
  str = StringUtil::convertToString(color_);
  xml_doc->addElement(element, "unexpected_change_color", str.c_str());
  xml_doc->addElement(element, "segment_type", segment_type_.c_str());
  if (ptm_ptr_ != nullptr) {
    ptm_ptr_->appendAbbrNameToXml(xml_doc, element);
  }
  xml_doc->addElement(element, "occurence", occu_.c_str());
  parent->appendChild(element);
}

}

