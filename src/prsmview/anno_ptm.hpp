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

#ifndef TOPPIC_PRSM_VIEW_ANNO_PTM_HPP_
#define TOPPIC_PRSM_VIEW_ANNO_PTM_HPP_

#include "common/xml/xml_dom_document.hpp"
#include "common/base/ptm.hpp"
#include "seq/mass_shift_type.hpp"
#include "prsmview/anno_ptm_position.hpp"

namespace toppic {

class AnnoPtm;
typedef std::shared_ptr<AnnoPtm> AnnoPtmPtr;
typedef std::vector<AnnoPtmPtr> AnnoPtmPtrVec;

class AnnoPtm {
 public:
  AnnoPtm(PtmPtr ptm_ptr, MassShiftTypePtr type_ptr);

  PtmPtr getPtmPtr() {return ptm_ptr_;}

  MassShiftTypePtr getTypePtr() {return type_ptr_;}

  void appendXml(XmlDOMDocument* xml_doc, xercesc::DOMElement* parent);

  void addOccurence(int left_pos, int right_pos, std::string anno);

  static AnnoPtmPtr findPtm(const AnnoPtmPtrVec &ptm_ptrs, PtmPtr ptm_ptr,
                            MassShiftTypePtr type_ptr);

 private:
  PtmPtr ptm_ptr_;

  MassShiftTypePtr type_ptr_;

  AnnoPtmPositionPtrVec occurences_;
};

}  // namespace toppic

#endif

