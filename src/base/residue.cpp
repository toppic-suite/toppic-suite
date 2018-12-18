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


#include <string>

#include <util/logger.hpp>

#include "base/amino_acid_base.hpp"
#include "base/ptm_base.hpp"
#include "base/residue.hpp"
#include "util/file_util.hpp"
#include "base/string_util.hpp"
#include "base/xml_dom_util.hpp"

namespace toppic {

Residue::Residue(xercesc::DOMElement* element) { 
  std::string acid_element_name = AminoAcid::getXmlElementName();
  xercesc::DOMElement* acid_element 
      = xml_dom_util::getChildElement(element, acid_element_name.c_str(), 0);
  acid_ptr_ = AminoAcidBase::getAminoAcidPtrFromXml(acid_element);
  std::string ptm_element_name = Ptm::getXmlElementName();
  xercesc::DOMElement* ptm_element 
      = xml_dom_util::getChildElement(element, ptm_element_name.c_str(), 0);
  ptm_ptr_ = PtmBase::getPtmPtrFromXml(ptm_element);  
  mass_ = acid_ptr_->getMonoMass() + ptm_ptr_->getMonoMass();
}

std::string Residue::toString(const std::string &delim_bgn, 
                              const std::string &delim_end) {
  if (PtmBase::isEmptyPtmPtr(ptm_ptr_)) {
    return acid_ptr_->getOneLetter();
  } else {
    return acid_ptr_->getOneLetter() + delim_bgn + ptm_ptr_->getAbbrName()
        + delim_end;
  }
}

void Residue::appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent,
                        const std::string & element_name) {
  xercesc::DOMElement* element = xml_doc->createElement(element_name.c_str());
  std::string str = string_util::convertToString(mass_);
  acid_ptr_->appendNameToXml(xml_doc, element);
  ptm_ptr_->appendAbbrNameToXml(xml_doc, element);
  parent->appendChild(element);
}

void Residue::appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent) {
  std::string element_name = Residue::getXmlElementName();
  appendXml(xml_doc, parent, element_name);
}

}  // namespace toppic
