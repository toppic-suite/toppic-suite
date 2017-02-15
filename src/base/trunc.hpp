// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#ifndef PROT_BASE_TRUNC_HPP_
#define PROT_BASE_TRUNC_HPP_

#include <string>
#include "base/residue.hpp"
#include "base/xml_dom_document.hpp"

namespace prot {

class Trunc {
 public:
  Trunc(const std::string &name, int trunc_len, 
        const std::string &trunc_residues,
        const std::string &allow_first_remain_residues_);

  Trunc(xercesc::DOMElement* element); 

  const std::string& getName() {return name_;}

  int getTruncLen() {return trunc_len_;}

  const ResiduePtrVec& getTruncResiduePtrVec() {return trunc_residue_ptr_vec_;}

  const ResiduePtrVec& getAllowFirstRemainResiduePtrs() {return allow_first_remain_residue_ptrs_;}

  double getShift() {return shift_;}

  static std::string getNameFromXml(xercesc::DOMElement * element);

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
