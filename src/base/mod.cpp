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
#include "base/mod.hpp"
#include "base/residue_base.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

Mod::Mod(ResiduePtr ori_residue_ptr, ResiduePtr mod_residue_ptr):
    ori_residue_ptr_(ori_residue_ptr),
    mod_residue_ptr_(mod_residue_ptr) {
    }

Mod::Mod(xercesc::DOMElement* element) { 
  xercesc::DOMElement* ori_residue_element 
      = XmlDomUtil::getChildElement(element, "ori_residue", 0);
  ori_residue_ptr_ = ResidueBase::getResiduePtrFromXml(ori_residue_element);
  xercesc::DOMElement* mod_residue_element 
      = XmlDomUtil::getChildElement(element, "mod_residue", 0);
  mod_residue_ptr_ = ResidueBase::getResiduePtrFromXml(mod_residue_element);
}

void Mod::appendToXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent){
  std::string element_name = Mod::getXmlElementName();
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  ori_residue_ptr_->appendXml(xml_doc, element, "ori_residue");
  mod_residue_ptr_->appendXml(xml_doc, element, "mod_residue");
  parent->appendChild(element);
}

}
