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


#include <string>

#include "common/util/logger.hpp"
#include "common/xml/xml_dom_document.hpp"
#include "common/xml/xml_dom_util.hpp"
#include "common/base/ptm_base.hpp"
#include "common/base/residue_base.hpp"
#include "common/base/n_term_mod.hpp"

namespace toppic {

NTermMod::NTermMod(ResiduePtr residue_ptr, PtmPtr ori_ptm_ptr, PtmPtr mod_ptm_ptr):
    residue_ptr_(residue_ptr),
    ori_ptm_ptr_(ori_ptm_ptr),
    mod_ptm_ptr_(mod_ptm_ptr) {}  

NTermMod::NTermMod(XmlDOMElement* element) {
  XmlDOMElement* residue_element
      = xml_dom_util::getChildElement(element, "residue", 0);
  residue_ptr_ = ResidueBase::getResiduePtrFromXml(residue_element);
  XmlDOMElement* ori_ptm_element
      = xml_dom_util::getChildElement(element, "ori_ptm", 0);
  ori_ptm_ptr_ = PtmBase::getPtmPtrFromXml(ori_ptm_element);
  XmlDOMElement* mod_ptm_element
      = xml_dom_util::getChildElement(element, "mod_ptm", 0);
  mod_ptm_ptr_ = PtmBase::getPtmPtrFromXml(mod_ptm_element);
}

bool NTermMod::isSame(NTermModPtr mod_ptr) {
  return residue_ptr_ == mod_ptr->getResiduePtr()
      && ori_ptm_ptr_ == mod_ptr->getOriPtmPtr()
      && mod_ptm_ptr_ == mod_ptr->getModPtmPtr();
}

double NTermMod::getReplaceShift() {
  if (PtmBase::isEmptyPtmPtr(ori_ptm_ptr_)) {
    return mod_ptm_ptr_->getMonoMass();
  }
  else {
    return mod_ptm_ptr_->getMonoMass() - ori_ptm_ptr_->getMonoMass();
  }
}

double NTermMod::getShift() {
  return mod_ptm_ptr_->getMonoMass();
}

void NTermMod::appendToXml(XmlDOMDocument* xml_doc, XmlDOMElement* parent) {
  std::string element_name = NTermMod::getXmlElementName();
  XmlDOMElement* element = xml_doc->createElement(element_name.c_str());
  residue_ptr_->appendXml(xml_doc, element, "residue");
  ori_ptm_ptr_->appendAbbrNameToXml(xml_doc, element, "ori_ptm");
  mod_ptm_ptr_->appendAbbrNameToXml(xml_doc, element, "mod_ptm");
  parent->appendChild(element);
}

}  // namespace toppic
