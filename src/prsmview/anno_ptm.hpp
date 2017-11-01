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


#ifndef PROT_ANNO_PTM_HPP_
#define PROT_ANNO_PTM_HPP_

#include "base/ptm.hpp"
#include "base/change_type.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

class AnnoPtm;
typedef std::shared_ptr<AnnoPtm> AnnoPtmPtr;
typedef std::vector<AnnoPtmPtr> AnnoPtmPtrVec;

class AnnoPtm {
 public:

  AnnoPtm(PtmPtr ptm_ptr, ChangeTypePtr change_type_ptr);

  PtmPtr getPtmPtr() {return ptm_ptr_;}

  ChangeTypePtr getChangeTypePtr() {return change_type_ptr_;}

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  void addOccurence(int pos, const std::string &acid_letter);

  static AnnoPtmPtr findPtm(const AnnoPtmPtrVec &ptm_ptrs, PtmPtr ptm_ptr, 
                            ChangeTypePtr change_type_ptr);

 private:
  PtmPtr ptm_ptr_;
  ChangeTypePtr change_type_ptr_;
  std::vector<std::pair<int, std::string>> occurences_;
};

}
#endif

