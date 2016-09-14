// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
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


#include "base/string_util.hpp"
#include "prsmview/anno_ptm.hpp"

namespace prot {

AnnoPtm::AnnoPtm(PtmPtr ptm_ptr, ChangeTypePtr change_type_ptr) {
  ptm_ptr_ = ptm_ptr;
  change_type_ptr_ = change_type_ptr;
}

void AnnoPtm::addOccurence(int pos, const std::string &acid_letter) {
  std::pair<int, std::string> new_occurence(pos, acid_letter);
  occurences_.push_back(new_occurence);
}

void AnnoPtm::appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent){
  xercesc::DOMElement* element = xml_doc->createElement("expected_change");
  std::string str = change_type_ptr_->getName();
  xml_doc->addElement(element, "change_type", str.c_str());
  ptm_ptr_->appendAbbrNameToXml(xml_doc,element);

  for (size_t i = 0; i < occurences_.size(); i++) {
    xercesc::DOMElement* position_element = xml_doc->createElement("occurence");
    std::string str = StringUtil::convertToString(occurences_[i].first);
    xml_doc->addElement(position_element, "position", str.c_str());
    xml_doc->addElement(position_element, "acid_letter", occurences_[i].second.c_str());
    element->appendChild(position_element);
  }
  parent->appendChild(element);
}

AnnoPtmPtr AnnoPtm::findPtm(const AnnoPtmPtrVec &ptm_ptrs, PtmPtr ptm_ptr, 
                            ChangeTypePtr change_type_ptr) {
  for (size_t i = 0; i < ptm_ptrs.size(); i++) {
    if ((ptm_ptrs[i]->getChangeTypePtr() == change_type_ptr) &&
        (ptm_ptrs[i]->getPtmPtr() == ptm_ptr)) {
      return ptm_ptrs[i];
    }
  }
  return nullptr;
}

}
