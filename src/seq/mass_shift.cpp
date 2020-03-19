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

#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "seq/mass_shift.hpp"

namespace toppic {

MassShift::MassShift(AlterPtr alter_ptr) {
  shift_ = alter_ptr->getMass();
  alter_vec_.push_back(alter_ptr);
  left_bp_pos_ = alter_vec_[0]->getLeftBpPos();
  right_bp_pos_ = alter_vec_[0]->getRightBpPos();
}

MassShift::MassShift(MassShiftPtr shift_ptr, int start) {
  shift_ = shift_ptr->getMassShift();
  left_bp_pos_ = shift_ptr->getLeftBpPos() - start;
  right_bp_pos_ = shift_ptr->getRightBpPos() - start;
  AlterPtrVec alter_ptrs = shift_ptr->getAlterPtrVec();
  for (size_t k = 0; k < alter_ptrs.size(); k++) {
    AlterPtr alter_ptr = Alter::geneAlterPtr(alter_ptrs[k], start); 
    alter_vec_.push_back(alter_ptr);
  }
}

MassShift::MassShift(int left_bp_pos, int right_bp_pos, 
                     double shift):
    left_bp_pos_(left_bp_pos),
    right_bp_pos_(right_bp_pos),
    shift_(shift) {}

AlterTypePtr MassShift::getTypePtr() {
  if (alter_vec_.size() == 0) {
    LOG_ERROR("Alter vector is empty!");
    exit(EXIT_FAILURE);
  } 
  return alter_vec_[0]->getTypePtr();
} 

MassShift::MassShift(XmlDOMElement* element) {
  left_bp_pos_ = xml_dom_util::getIntChildValue(element, "left_bp_pos", 0);
  right_bp_pos_ = xml_dom_util::getIntChildValue(element, "right_bp_pos", 0);

  shift_ = xml_dom_util::getDoubleChildValue(element, "shift", 0);

  std::string alter_element_name = Alter::getXmlElementName();
  std::string alter_element_list = alter_element_name + "_list";
  XmlDOMElement* alter_list_element 
      = xml_dom_util::getChildElement(element, alter_element_list.c_str(), 0);

  int alter_len = xml_dom_util::getChildCount(alter_list_element, 
                                              alter_element_name.c_str());
  for (int i = 0; i < alter_len; i++) {
    XmlDOMElement* alter_element
        = xml_dom_util::getChildElement(alter_list_element, alter_element_name.c_str(), i);
    alter_vec_.push_back(std::make_shared<Alter>(alter_element));
  }
}

std::string MassShift::getAnnoStr() {
  std::string seq_str;
  if (getTypePtr() == AlterType::UNEXPECTED) {
    if (alter_vec_[0]->getLocalAnno() != nullptr) {
      seq_str = alter_vec_[0]->getLocalAnno()->getPtmPtr()->getAbbrName();
    } else {
      seq_str = str_util::fixedToString(shift_, 4);
    }
  } else {
    for (size_t i = 0; i < alter_vec_.size(); i++) {
      seq_str += alter_vec_[i]->getModPtr()->getModResiduePtr()->getPtmPtr()->getAbbrName();
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
  xml_doc->addElement(element, "left_bp_pos", str.c_str());
  str = str_util::toString(right_bp_pos_);
  xml_doc->addElement(element, "right_bp_pos", str.c_str());
  str = str_util::toString(shift_);
  xml_doc->addElement(element, "shift", str.c_str());

  element_name = Alter::getXmlElementName() + "_list";
  XmlDOMElement* cl = xml_doc->createElement(element_name.c_str());
  for (size_t i = 0; i < alter_vec_.size(); i++) {
    alter_vec_[i]->appendXml(xml_doc, cl);
  }
  element->appendChild(cl);
  parent->appendChild(element);
}

bool MassShift::cmpPosInc(const MassShiftPtr &a, const MassShiftPtr &b) {
  if (a->getLeftBpPos() < b->getLeftBpPos()) {
    return true;
  } else if (a->getLeftBpPos() > b->getLeftBpPos()) {
    return false;
  } else {
    return a->getRightBpPos() < b->getRightBpPos();
  }
}

}  // namespace toppic
