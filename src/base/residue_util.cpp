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
#include "base/acid_base.hpp"
#include "base/ptm_base.hpp"
#include "base/fasta_seq.hpp"
#include "base/residue_base.hpp"
#include "base/residue_util.hpp"

namespace prot {

ResiduePtrVec ResidueUtil::convertStrToResiduePtrVec(const std::string &seq) {
  ResiduePtrVec residue_ptr_vec;
  std::string seq2 = FastaSeq::rmChar(seq);
  for (size_t i = 0; i < seq2.length(); i++) {
    AcidPtr acid_ptr = AcidBase::getAcidPtrByOneLetter(seq2.substr(i, 1));
    PtmPtr ptm_ptr = PtmBase::getEmptyPtmPtr();
    ResiduePtr residue_ptr = ResidueBase::getBaseResiduePtr(acid_ptr, ptm_ptr);
    residue_ptr_vec.push_back(residue_ptr);
  }
  return residue_ptr_vec;
}

void applyFixedMod(ResiduePtrVec &residue_ptrs, const ModPtrVec &fix_mod_ptr_vec) {
  for (size_t i = 0; i < residue_ptrs.size(); i++) {
    for (size_t j = 0; j < fix_mod_ptr_vec.size(); j++) {
      if (residue_ptrs[i] == fix_mod_ptr_vec[j]->getOriResiduePtr()) {
        residue_ptrs[i] = fix_mod_ptr_vec[j]->getModResiduePtr();
        break;
      }
    }
  }
}

ResiduePtrVec ResidueUtil::convertStrToResiduePtrVec(const std::string &seq, 
                                                     const ModPtrVec &fix_mod_ptr_vec) {
  ResiduePtrVec residue_ptrs = ResidueUtil::convertStrToResiduePtrVec(seq);
  applyFixedMod(residue_ptrs, fix_mod_ptr_vec);
  return residue_ptrs;
}

ResiduePtrVec ResidueUtil::convertStrToResiduePtrVec(const StringPairVec &string_pair_vec) {
  ResiduePtrVec residue_ptr_vec;
  for (size_t i = 0; i < string_pair_vec.size(); i++) {
    std::string acid_str = string_pair_vec[i].first;
    AcidPtr acid_ptr = AcidBase::getAcidPtrByOneLetter(acid_str);
    std::string ptm_str = string_pair_vec[i].second;
    PtmPtr ptm_ptr = PtmBase::getPtmPtrByAbbrName(ptm_str);
    ResiduePtr residue_ptr = ResidueBase::getBaseResiduePtr(acid_ptr, ptm_ptr);
    residue_ptr_vec.push_back(residue_ptr);
  }
  return residue_ptr_vec;
}

ResiduePtrVec ResidueUtil::convertStrToResiduePtrVec(const StringPairVec &string_pair_vec,  
                                                     const ModPtrVec &fix_mod_ptr_vec) {
  ResiduePtrVec residue_ptrs = ResidueUtil::convertStrToResiduePtrVec(string_pair_vec);
  applyFixedMod(residue_ptrs, fix_mod_ptr_vec);
  return residue_ptrs;
}

int ResidueUtil::findResidue(const ResiduePtrVec &residue_list, 
                             ResiduePtr residue_ptr) {
  for (size_t i = 0; i < residue_list.size(); i++) {
    if (residue_list[i] == residue_ptr) {
      return i;
    }
  }
  return -1;
}

double ResidueUtil::compResiduePtrVecMass(const ResiduePtrVec &ptr_vec) {
  double mass = 0;
  for (size_t i = 0; i < ptr_vec.size(); i++) {
    mass += ptr_vec[i]->getMass();
  }
  return mass;
}

}
