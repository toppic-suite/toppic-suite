// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "seq/alter.hpp"

#include <memory>
#include <string>

#include "common/base/mod_base.hpp"
#include "common/util/str_util.hpp"
#include "common/xml/xml_dom_document.hpp"
#include "common/xml/xml_dom_util.hpp"

namespace toppic {

Alter::Alter(int left_bp_pos, int right_bp_pos, const AlterTypePtr& type_ptr,
             double mass, const ModPtr& mod_ptr)
    : left_bp_pos_(left_bp_pos),
      right_bp_pos_(right_bp_pos),
      type_ptr_(type_ptr),
      mass_(mass),
      mod_ptr_(mod_ptr),
      local_anno_ptr_(nullptr) {}

Alter::Alter(XmlDOMElement element) {
  left_bp_pos_ = xml_dom_util::getIntChildValue(element, "left_bp_pos", 0);
  right_bp_pos_ = xml_dom_util::getIntChildValue(element, "right_bp_pos", 0);
  std::string type_element_name = AlterType::getXmlElementName();
  XmlDOMElement type_element =
      xml_dom_util::getChildElement(element, type_element_name.c_str(), 0);
  type_ptr_ = AlterType::getTypePtrFromXml(type_element);
  mass_ = xml_dom_util::getDoubleChildValue(element, "mass", 0);
  std::string mod_element_name = Mod::getXmlElementName();

  int mod_count =
      xml_dom_util::getChildCount(element, mod_element_name.c_str());
  if (mod_count != 0) {
    XmlDOMElement mod_element =
        xml_dom_util::getChildElement(element, mod_element_name.c_str(), 0);
    mod_ptr_ = ModBase::getModPtrFromXml(mod_element);
  }

  std::string local_element_name = LocalAnno::getXmlElementName();
  int local_count =
      xml_dom_util::getChildCount(element, local_element_name.c_str());
  if (local_count != 0) {
    XmlDOMElement local_element =
        xml_dom_util::getChildElement(element, local_element_name.c_str(), 0);
    local_anno_ptr_ = std::make_shared<LocalAnno>(local_element);
  }
}

void Alter::appendXml(XmlDOMDocument* xml_doc, XmlDOMElement parent) const {
  std::string element_name = Alter::getXmlElementName();
  XmlDOMElement element = xml_doc->addElement(parent, element_name.c_str());
  std::string str = std::to_string(left_bp_pos_);
  xml_doc->addElement(element, "left_bp_pos", str.c_str());
  str = std::to_string(right_bp_pos_);
  xml_doc->addElement(element, "right_bp_pos", str.c_str());
  type_ptr_->appendXml(xml_doc, element);
  str = str_util::toString(mass_);
  xml_doc->addElement(element, "mass", str.c_str());
  if (mod_ptr_ != nullptr) {
    mod_ptr_->appendToXml(xml_doc, element);
  }
  if (local_anno_ptr_ != nullptr) {
    local_anno_ptr_->appendToXml(xml_doc, element);
  }
}

AlterPtr Alter::genAlterPtr(const AlterPtr& ori_ptr, int start_pos) {
  int left_bp_pos = ori_ptr->left_bp_pos_ - start_pos;
  int right_bp_pos = ori_ptr->right_bp_pos_ - start_pos;
  AlterTypePtr type_ptr = ori_ptr->type_ptr_;
  double mass = ori_ptr->mass_;
  ModPtr mod_ptr = ori_ptr->mod_ptr_;
  AlterPtr alter_ptr = std::make_shared<Alter>(left_bp_pos, right_bp_pos,
                                               type_ptr, mass, mod_ptr);
  alter_ptr->local_anno_ptr_ = ori_ptr->local_anno_ptr_;
  return alter_ptr;
}

void Alter::setLocalAnno(const LocalAnnoPtr& p) {
  local_anno_ptr_ = p;
  if (p != nullptr) {
    left_bp_pos_ = p->getLeftBpPos();
    right_bp_pos_ = p->getRightBpPos() + 1;
  }
}

}  // namespace toppic
