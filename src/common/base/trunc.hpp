//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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

#ifndef TOPPIC_COMMON_BASE_TRUNC_HPP_
#define TOPPIC_COMMON_BASE_TRUNC_HPP_

#include "common/xml/xml_dom_element.hpp"

#include "common/base/residue.hpp"

namespace toppic {

class Trunc {
 public:
  Trunc(const std::string &name, int trunc_len, 
        const std::string &trunc_residues,
        const std::string &allow_first_remain_residues_);

  Trunc(XmlDOMElement* element); 

  const std::string& getName() {return name_;}

  int getTruncLen() {return trunc_len_;}

  const ResiduePtrVec& getTruncResiduePtrVec() {return trunc_residue_ptr_vec_;}

  const ResiduePtrVec& getAllowFirstRemainResiduePtrs() {return allow_first_remain_residue_ptrs_;}

  double getShift() {return shift_;}

  static std::string getNameFromXml(XmlDOMElement * element);

  static std::string getXmlElementName() {return "truncation";}

 private:
  std::string name_;
  int trunc_len_;
  ResiduePtrVec trunc_residue_ptr_vec_;
  ResiduePtrVec allow_first_remain_residue_ptrs_;
  double shift_;
};

typedef std::shared_ptr<Trunc> TruncPtr;
typedef std::vector<TruncPtr> TruncPtrVec;

}

#endif
