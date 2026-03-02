//Copyright (c) 2014 - 2025, The Trustees of Indiana University, Tulane University.
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

#include <stdexcept>
#include <string>

#include "common/util/logger.hpp"
#include "common/xml/xml_dom_document.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "common/base/ptm_base.hpp"
#include "common/base/mod.hpp"
#include "common/base/residue_base.hpp"

namespace toppic {

Mod::Mod(ResiduePtr ori_residue_ptr, ResiduePtr mod_residue_ptr, ModTypePtr mod_type_ptr):
    ori_residue_ptr_(ori_residue_ptr),
    mod_residue_ptr_(mod_residue_ptr),
    mod_type_ptr_(mod_type_ptr) {
}

Mod::Mod(XmlDOMElement* element) {
  XmlDOMElement* ori_residue_element
      = xml_dom_util::getChildElement(element, "ori_residue", 0);
  ori_residue_ptr_ = ResidueBase::getResiduePtrFromXml(ori_residue_element);
  XmlDOMElement* mod_residue_element
      = xml_dom_util::getChildElement(element, "mod_residue", 0);
  mod_residue_ptr_ = ResidueBase::getResiduePtrFromXml(mod_residue_element);
  std::string mod_type_name
    = xml_dom_util::getChildValue(element, "mod_type", 0);
  mod_type_ptr_ = ModType::getModTypePtrByName(mod_type_name);
  if (mod_type_ptr_ == nullptr) {
    throw std::runtime_error("Unknown mod type: " + mod_type_name);
  }
}

bool Mod::isSame(ModPtr mod_ptr) const {
  return ori_residue_ptr_ == mod_ptr->getOriResiduePtr()
      && mod_residue_ptr_ == mod_ptr->getModResiduePtr()
      && mod_type_ptr_ == mod_ptr->getModTypePtr();
}

double Mod::getReplaceShift() const {
  if (PtmBase::isEmptyPtmPtr(ori_residue_ptr_->getPtmPtr())) {
    return mod_residue_ptr_->getPtmPtr()->getMonoMass();
  }
  else {
    return mod_residue_ptr_->getMass() - ori_residue_ptr_->getMass();
  }
}

double Mod::getShift() const {
  return mod_residue_ptr_->getPtmPtr()->getMonoMass(); 
}

void Mod::appendToXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent) const {
  std::string element_name = Mod::getXmlElementName();
  XmlDOMElement* element = xml_doc->createElement(element_name.c_str());
  ori_residue_ptr_->appendXml(xml_doc, element, "ori_residue");
  mod_residue_ptr_->appendXml(xml_doc, element, "mod_residue");
  xml_doc->addElement(element, "mod_type", mod_type_ptr_->getName().c_str());
  parent->appendChild(element);
}

}  // namespace toppic
