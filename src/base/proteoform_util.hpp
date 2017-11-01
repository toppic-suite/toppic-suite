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
