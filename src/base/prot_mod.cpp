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


#include "base/logger.hpp"
#include "base/ptm_base.hpp"
#include "base/trunc_base.hpp"
#include "base/mod_base.hpp"
#include "base/prot_mod.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

ProtMod::ProtMod(const std::string &name, const std::string &type,
                 TruncPtr trunc_ptr, ModPtr mod_ptr): 
    name_(name),
    type_(type),
    trunc_ptr_(trunc_ptr),
    mod_ptr_(mod_ptr) {
      mod_pos_ = trunc_ptr->getTruncLen();
      prot_shift_ = trunc_ptr_->getShift() + mod_ptr_->getShift();
      pep_shift_ = mod_ptr_->getShift();
    }

ProtMod::ProtMod(xercesc::DOMElement* element) { 
  name_ = xml_dom_util::getChildValue(element, "name", 0);
  type_ = xml_dom_util::getChildValue(element, "type", 0);
  std::string trunc_element_name = Trunc::getXmlElementName();
  xercesc::DOMElement* trunc_element 
      = xml_dom_util::getChildElement(element, trunc_element_name.c_str(), 0);
  trunc_ptr_ = TruncBase::getTruncPtrFromXml(trunc_element);
  std::string mod_element_name = Mod::getXmlElementName();
  xercesc::DOMElement* mod_element 
      = xml_dom_util::getChildElement(element, mod_element_name.c_str(), 0);
  mod_ptr_= ModBase::getModPtrFromXml(mod_element); 
  mod_pos_ = trunc_ptr_->getTruncLen();
  prot_shift_ = trunc_ptr_->getShift() + mod_ptr_->getShift();
  pep_shift_ = mod_ptr_->getShift();
}

void ProtMod::appendNameToXml(XmlDOMDocument* xml_doc,
                              xercesc::DOMElement* parent){
  std::string element_name = ProtMod::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  xml_doc->addElement(element, "name", name_.c_str());
  parent->appendChild(element);
}

std::string ProtMod::getNameFromXml(xercesc::DOMElement * element) {
  std::string name = xml_dom_util::getChildValue(element, "name", 0);
  return name;
}

bool ProtMod::isAcetylation() {
  if (mod_ptr_->getModResiduePtr()->getPtmPtr() == PtmBase::getPtmPtr_Acetylation()) {
    return true;
  }
  else {
    return false;
  }
}

}
