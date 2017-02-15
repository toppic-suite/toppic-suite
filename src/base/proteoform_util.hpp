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


#ifndef PROT_BASE_PROTEOFORM_UTIL_HPP_
#define PROT_BASE_PROTEOFORM_UTIL_HPP_

#include "base/proteoform.hpp"


namespace prot {

class ProteoformUtil {
 public:
  /* calculate frequencies for n_terminal_residues */
  static ResFreqPtrVec compNTermResidueFreq(const ProteoformPtrVec &prot_mod_forms);

  /* calculater frequences for all residues */
  static ResFreqPtrVec compResidueFreq(const ResiduePtrVec &residue_list,
                                       const ProteoformPtrVec &raw_mods);

  static bool isSameSeqAndMass(ProteoformPtr a, ProteoformPtr b, double ppo);

  static bool isStrictCompatiablePtmSpecies(ProteoformPtr a, ProteoformPtr b, double ppo);

  static ProteoformPtrVec2D divideProteoIntoBlocks(const ProteoformPtrVec &proteo_ptrs, 
                                                   int db_block_size);

  static std::vector<double> getNTermShift(ProteoformPtr db_form_ptr,
                                        const ProtModPtrVec &prot_mod_ptrs);

  static std::vector<std::vector<double>> getNTermShift2D(
      ProteoformPtrVec db_form_ptr_vec, const ProtModPtrVec &prot_mod_ptrs);

  static std::vector<double> getNTermAcets(ProteoformPtr db_form_ptr,
                                        const ProtModPtrVec &prot_mod_ptrs);

  static std::vector<std::vector<double>> getNTermAcet2D(
      ProteoformPtrVec db_form_ptr_vec, const ProtModPtrVec &prot_mod_ptrs);

};

} /* namespace prot */

#endif 
