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

#include <cmath>
#include <utility>
#include <string>

#include "util/string_util.hpp"
#include "prsmview/anno_segment.hpp"

namespace toppic {

void AnnoSegment::addOccurence(int pos, const std::string &acid_letter) {
  std::pair<int, std::string> new_occurence(pos + 1, acid_letter);
  occurences_.push_back(new_occurence);
}

std::string AnnoSegment::getResidueAnno() {
  if (anno_ != "") {
    return anno_;
  }

  if (ptm_ptr_ != nullptr) {
    anno_ += "PTM: " + ptm_ptr_->getName() + "\n";
    for (size_t i = 0; i < occurences_.size(); i++) {
      if (score_[i] == 1) score_[i] = 0.999;

      if (score_[i] == 0) continue;

      score_[i] = std::floor(score_[i] * 1000) / 10;

      anno_ += "Site: " + occurences_[i].second + string_util::convertToString(occurences_[i].first) + " ";
      anno_ += "Confidence: " + string_util::convertToString(score_[i], 1) + "%\n";
    }
  }
  return anno_;
}

void AnnoSegment::appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent) {
  xercesc::DOMElement* element;
  if (mass_shift_type_ == toppic::MassShiftType::UNEXPECTED) {
    if (ptm_ptr_ != nullptr) {
      element = xml_doc->createElement("characterized_change");
    } else {
      element = xml_doc->createElement("unexpected_change");
    }
    std::string str = string_util::convertToString(left_pos_);
    xml_doc->addElement(element, "left_position", str.c_str());

    str = string_util::convertToString(right_pos_);
    xml_doc->addElement(element, "right_position", str.c_str());

    xml_doc->addElement(element, "match_seq", match_seq_.c_str());

    str = string_util::convertToString(color_);
    xml_doc->addElement(element, "unexpected_change_color", str.c_str());

    xml_doc->addElement(element, "segment_type", segment_type_.c_str());

    xml_doc->addElement(element, "anno", getResidueAnno().c_str());

    if (ptm_ptr_ != nullptr) {
      ptm_ptr_->appendAbbrNameToXml(xml_doc, element);
    }

    std::string occu;

    if (occurences_.size() == 1) {
      occu = occurences_[0].second + string_util::convertToString(occurences_[0].first);
    } else if (occurences_.size() > 1) {
      occu = occurences_[0].second + string_util::convertToString(occurences_[0].first);
      occu += " - ";
      occu += occurences_[occurences_.size() - 1].second
          + string_util::convertToString(occurences_[occurences_.size() - 1].first);
    }

    xml_doc->addElement(element, "occurence", occu.c_str());
  } else {
    element = xml_doc->createElement("variable_change");
    std::string str = string_util::convertToString(left_pos_);
    xml_doc->addElement(element, "left_position", str.c_str());

    str = string_util::convertToString(right_pos_);
    xml_doc->addElement(element, "right_position", str.c_str());

    xml_doc->addElement(element, "match_seq", match_seq_.c_str());

    str = string_util::convertToString(color_);
    xml_doc->addElement(element, "unexpected_change_color", str.c_str());

    xml_doc->addElement(element, "segment_type", segment_type_.c_str());

    std::string occu;

    if (occurences_.size() == 1) {
      occu = occurences_[0].second + string_util::convertToString(occurences_[0].first);
    } else if (occurences_.size() > 1) {
      occu = occurences_[0].second + string_util::convertToString(occurences_[0].first);
      occu += " - ";
      occu += occurences_[occurences_.size() - 1].second
          + string_util::convertToString(occurences_[occurences_.size() - 1].first);
    }

    xml_doc->addElement(element, "occurence", occu.c_str());

    xml_doc->addElement(element, "anno", getResidueAnno().c_str());
  }
  parent->appendChild(element);
}

}  // namespace toppic

