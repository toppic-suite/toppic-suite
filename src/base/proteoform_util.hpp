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


#ifndef PROT_BASE_PROTEOFORM_UTIL_HPP_
#define PROT_BASE_PROTEOFORM_UTIL_HPP_

#include <vector>

#include "base/proteoform.hpp"

namespace prot {

namespace proteoform_util {

// calculate frequencies for n_terminal_residues
ResFreqPtrVec compNTermResidueFreq(const ProteoformPtrVec &prot_mod_forms);

// calculater frequences for all residues
ResFreqPtrVec compResidueFreq(const ResiduePtrVec &residue_list,
                              const ProteoformPtrVec &raw_mods);

bool isSameSeqAndMass(ProteoformPtr a, ProteoformPtr b, double ppo);

bool isStrictCompatiablePtmSpecies(ProteoformPtr a, ProteoformPtr b, double ppo);

ProteoformPtrVec2D divideProteoIntoBlocks(const ProteoformPtrVec &proteo_ptrs,
                                          int db_block_size);

std::vector<double> getNTermShift(ProteoformPtr db_form_ptr,
                                  const ProtModPtrVec &prot_mod_ptrs);

std::vector<std::vector<double> > getNTermShift2D(const ProteoformPtrVec & db_form_ptr_vec,
                                                  const ProtModPtrVec & prot_mod_ptrs);

std::vector<double> getNTermAcets(ProteoformPtr db_form_ptr,
                                  const ProtModPtrVec & prot_mod_ptrs);

std::vector<std::vector<double> > getNTermAcet2D(const ProteoformPtrVec & db_form_ptr_vec,
                                                 const ProtModPtrVec & prot_mod_ptrs);

}  // namespace proteoform_util

}  // namespace prot

#endif
