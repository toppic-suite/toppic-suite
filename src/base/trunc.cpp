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


#include "base/logger.hpp"
#include "base/trunc.hpp"
#include "base/residue_util.hpp"
#include "base/xml_dom_util.hpp"

namespace prot {

Trunc::Trunc(const std::string &name, int trunc_len, 
             const std::string &trunc_residues,
             const std::string &allow_first_remain_residues) {
  name_ = name;
  trunc_len_ = trunc_len;
  trunc_residue_ptr_vec_ = ResidueUtil::convertStrToResiduePtrVec(trunc_residues);
  allow_first_remain_residue_ptrs_ = ResidueUtil::convertStrToResiduePtrVec(allow_first_remain_residues);
  shift_ = -ResidueUtil::compResiduePtrVecMass(trunc_residue_ptr_vec_);
}

Trunc::Trunc(xercesc::DOMElement* element) { 
  name_ = XmlDomUtil::getChildValue(element, "name", 0);
  trunc_len_ = XmlDomUtil::getIntChildValue(element, "trunc_len", 0);
  std::string trunc_residues = XmlDomUtil::getChildValue(element, "trunc_residues", 0);
  LOG_DEBUG( "name " << name_ << " str " << trunc_residues << " trunc len " << trunc_len_);
  trunc_residue_ptr_vec_ = ResidueUtil::convertStrToResiduePtrVec(trunc_residues);
  std::string allow_first_remain_residues = XmlDomUtil::getChildValue(element, "allow_first_remain_residues", 0);
  allow_first_remain_residue_ptrs_ = ResidueUtil::convertStrToResiduePtrVec(allow_first_remain_residues);
  shift_ = -ResidueUtil::compResiduePtrVecMass(trunc_residue_ptr_vec_);
}

std::string Trunc::getNameFromXml(xercesc::DOMElement * element) {
  std::string name = XmlDomUtil::getChildValue(element, "name", 0);
  return name;
}

}
