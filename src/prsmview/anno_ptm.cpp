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


#include <utility>
#include <string>

#include "common/util/str_util.hpp"
#include "prsmview/anno_ptm.hpp"

namespace toppic {

AnnoPtm::AnnoPtm(PtmPtr ptm_ptr, MassShiftTypePtr type_ptr) {
  ptm_ptr_ = ptm_ptr;
  type_ptr_ = type_ptr;
}

void AnnoPtm::addOccurence(int left_pos, int right_pos, const std::string anno) {
  AnnoPtmPositionPtr pos_ptr 
      = std::make_shared<AnnoPtmPosition>(left_pos, right_pos, anno);
  occurences_.push_back(pos_ptr);
}

void AnnoPtm::appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent) {
  xercesc::DOMElement* element = xml_doc->createElement("ptm");
  std::string str = type_ptr_->getName();
  xml_doc->addElement(element, "ptm_type", str.c_str());
  ptm_ptr_->appendAbbrNameToXml(xml_doc, element);

  for (size_t i = 0; i < occurences_.size(); i++) {
    occurences_[i]->appendXml(xml_doc, element);
  }
  parent->appendChild(element);
}

AnnoPtmPtr AnnoPtm::findPtm(const AnnoPtmPtrVec &ptm_ptrs, PtmPtr ptm_ptr,
                            MassShiftTypePtr type_ptr) {
  for (size_t i = 0; i < ptm_ptrs.size(); i++) {
    if ((ptm_ptrs[i]->getTypePtr() == type_ptr) &&
        (ptm_ptrs[i]->getPtmPtr() == ptm_ptr)) {
      return ptm_ptrs[i];
    }
  }
  return nullptr;
}

}  // namespace toppic
