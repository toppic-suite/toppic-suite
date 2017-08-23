// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <string>

#include "base/mod_base.hpp"
#include "base/change.hpp"
#include "base/string_util.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

Change::Change(xercesc::DOMElement* element) {
  left_bp_pos_ = XmlDomUtil::getIntChildValue(element, "left_bp_pos", 0);
  right_bp_pos_ = XmlDomUtil::getIntChildValue(element, "right_bp_pos", 0);
  std::string ct_element_name = ChangeType::getXmlElementName();
  xercesc::DOMElement* ct_element
      = XmlDomUtil::getChildElement(element, ct_element_name.c_str(), 0);
  change_type_ptr_ = ChangeType::getChangeTypePtrFromXml(ct_element);
  mass_shift_ = XmlDomUtil::getDoubleChildValue(element, "mass_shift", 0);
  std::string mod_element_name = Mod::getXmlElementName();
  int mod_count = XmlDomUtil::getChildCount(element, mod_element_name.c_str());
  if (mod_count != 0) {
    xercesc::DOMElement* mod_element
        = XmlDomUtil::getChildElement(element, mod_element_name.c_str(), 0);
    mod_ptr_ = ModBase::getModPtrFromXml(mod_element);
  }
  std::string local_element_name = LocalAnno::getXmlElementName();;
  int local_count = XmlDomUtil::getChildCount(element, local_element_name.c_str());
  if (local_count != 0) {
    xercesc::DOMElement * local_element
        = XmlDomUtil::getChildElement(element, local_element_name.c_str(), 0);
    local_anno_ptr_ = LocalAnnoPtr(new LocalAnno(local_element));
  }
}

void Change::appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent) {
  std::string element_name = Change::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  std::string str = StringUtil::convertToString(left_bp_pos_);
  xml_doc->addElement(element, "left_bp_pos", str.c_str());
  str = StringUtil::convertToString(right_bp_pos_);
  xml_doc->addElement(element, "right_bp_pos", str.c_str());
  change_type_ptr_->appendXml(xml_doc, element);
  str = StringUtil::convertToString(mass_shift_);
  xml_doc->addElement(element, "mass_shift", str.c_str());
  if (mod_ptr_ != nullptr) {
    mod_ptr_->appendToXml(xml_doc, element);
  }
  if (local_anno_ptr_ != nullptr) {
    local_anno_ptr_->appendToXml(xml_doc, element);
  }
  parent->appendChild(element);
}

bool Change::cmpPosInc(const ChangePtr &a, const ChangePtr &b) {
  if (a->getLeftBpPos() < b->getLeftBpPos()) {
    return true;
  } else if (a->getLeftBpPos() > b->getLeftBpPos()) {
    return false;
  } else {
    return a->getRightBpPos() < b->getRightBpPos();
  }
}

ChangePtr Change::geneChangePtr(ChangePtr ori_ptr, int start_pos) {
  int left_bp_pos = ori_ptr->left_bp_pos_ - start_pos;
  int right_bp_pos = ori_ptr->right_bp_pos_ - start_pos;
  ChangeTypePtr change_type_ptr = ori_ptr->change_type_ptr_;
  double mass_shift = ori_ptr->mass_shift_;
  ModPtr mod_ptr = ori_ptr->mod_ptr_;
  ChangePtr change_ptr(
      new Change(left_bp_pos, right_bp_pos, change_type_ptr, mass_shift, mod_ptr));
  return change_ptr;
}

void Change::setLocalAnno(LocalAnnoPtr p) {
  local_anno_ptr_ = p;
  if (p != nullptr) {
    left_bp_pos_ = p->getLeftBpPos();
    right_bp_pos_ = p->getRightBpPos() + 1;
  }
}

}  // namespace prot
