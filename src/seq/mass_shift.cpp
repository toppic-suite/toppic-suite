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

#include <string>

#include "common/util/str_util.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "seq/mass_shift.hpp"

namespace toppic {

MassShift::MassShift(int left_bp_pos, int right_bp_pos, MassShiftTypePtr type_ptr):
    left_bp_pos_(left_bp_pos),
    right_bp_pos_(right_bp_pos),
    type_ptr_(type_ptr),
    shift_(0.0) { }

MassShift::MassShift(XmlDOMElement* element) {
  left_bp_pos_ = xml_dom_util::getIntChildValue(element, "shift_left_bp_pos", 0);

  right_bp_pos_ = xml_dom_util::getIntChildValue(element, "shift_right_bp_pos", 0);

  std::string ct_element_name = MassShiftType::getXmlElementName();
  XmlDOMElement* ct_element
      = xml_dom_util::getChildElement(element, ct_element_name.c_str(), 0);
  type_ptr_ = MassShiftType::getChangeTypePtrFromXml(ct_element);

  shift_ = xml_dom_util::getDoubleChildValue(element, "shift", 0);

  XmlDOMElement* change_list_element = xml_dom_util::getChildElement(element, "change_list", 0);

  std::string change_element_name = Change::getXmlElementName();
  int change_len = xml_dom_util::getChildCount(change_list_element, change_element_name.c_str());

  for (int i = 0; i < change_len; i++) {
    XmlDOMElement* change_element
        = xml_dom_util::getChildElement(change_list_element, change_element_name.c_str(), i);
    change_vec_.push_back(std::make_shared<Change>(change_element));
  }
}

void MassShift::setChangePtr(ChangePtr change) {
  shift_ += change->getMass();
  change_vec_.push_back(change);
  left_bp_pos_ = change_vec_[0]->getLeftBpPos();
  right_bp_pos_ = change_vec_[0]->getRightBpPos();
  for (size_t k = 0; k < change_vec_.size(); k++) {
    if (change_vec_[k]->getLeftBpPos() < left_bp_pos_) {
      left_bp_pos_ = change_vec_[k]->getLeftBpPos();
    }
    if (change_vec_[k]->getRightBpPos() > right_bp_pos_) {
      right_bp_pos_ = change_vec_[k]->getRightBpPos();
    }
  }
}

std::string MassShift::getAnnoStr() {
  std::string seq_str;

  if (getTypePtr() == MassShiftType::UNEXPECTED) {
    if (change_vec_[0]->getLocalAnno() != nullptr) {
      seq_str = change_vec_[0]->getLocalAnno()->getPtmPtr()->getAbbrName();
    } else {
      seq_str = str_util::toString(shift_, 4);
    }
  } else {
    for (size_t i = 0; i < change_vec_.size(); i++) {
      seq_str += change_vec_[i]->getModPtr()->getModResiduePtr()->getPtmPtr()->getAbbrName();
      seq_str += ";";
    }
    seq_str.pop_back();
  }

  return seq_str;
}

void MassShift::appendXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent) {
  std::string element_name = getXmlElementName();
  XmlDOMElement* element = xml_doc->createElement(element_name.c_str());
  std::string str = str_util::toString(left_bp_pos_);
  xml_doc->addElement(element, "shift_left_bp_pos", str.c_str());
  str = str_util::toString(right_bp_pos_);
  xml_doc->addElement(element, "shift_right_bp_pos", str.c_str());
  type_ptr_->appendXml(xml_doc, element);
  str = str_util::toString(shift_);
  xml_doc->addElement(element, "shift", str.c_str());

  element_name = Change::getXmlElementName() + "_list";
  XmlDOMElement* cl = xml_doc->createElement(element_name.c_str());
  for (size_t i = 0; i < change_vec_.size(); i++) {
    change_vec_[i]->appendXml(xml_doc, cl);
  }
  element->appendChild(cl);
  parent->appendChild(element);
}

bool MassShift::cmpPosInc(const MassShiftPtr & a, const MassShiftPtr & b) {
  if (a->getLeftBpPos() < b->getLeftBpPos()) {
    return true;
  } else if (a->getLeftBpPos() > b->getLeftBpPos()) {
    return false;
  } else {
    return a->getRightBpPos() < b->getRightBpPos();
  }
}

}  // namespace toppic
