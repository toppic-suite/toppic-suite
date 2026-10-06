// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef TOPPIC_SEQ_PROTEOFORM_UTIL_HPP_
#define TOPPIC_SEQ_PROTEOFORM_UTIL_HPP_

#include <string>
#include <vector>

#include "seq/proteoform.hpp"

namespace toppic {

namespace proteoform_util {

// calculate frequencies for n_terminal_residues
ResFreqPtrVec compNTermResidueFreq(const ProteoformPtrVec& prot_mod_forms);

// calculate frequencies for all residues
ResFreqPtrVec compResidueFreq(const ResiduePtrVec& residue_list,
                              const ProteoformPtrVec& raw_mods);

bool isSameSeqAndMass(const ProteoformPtr& a, const ProteoformPtr& b,
                      double ppo);

bool isStrictCompatiablePtmSpecies(const ProteoformPtr& a,
                                   const ProteoformPtr& b, double ppo);

ProteoformPtrVec2D divideProteoIntoBlocks(const ProteoformPtrVec& proteo_ptrs,
                                          int db_block_size);

std::vector<double> getNTermShift(const ProteoformPtr& db_form_ptr,
                                  const ProtModPtrVec& prot_mod_ptrs);

std::vector<std::vector<double> > getNTermShift2D(
    const ProteoformPtrVec& db_form_ptr_vec,
    const ProtModPtrVec& prot_mod_ptrs);

std::vector<double> getNTermAcets(const ProteoformPtr& db_form_ptr,
                                  const ProtModPtrVec& prot_mod_ptrs);

std::vector<std::vector<double> > getNTermAcet2D(
    const ProteoformPtrVec& db_form_ptr_vec,
    const ProtModPtrVec& prot_mod_ptrs);

// The ProForma string of a proteoform (HUPO-PSI ProForma 2.1,
// https://github.com/HUPO-PSI/ProForma): its residues (first to last
// residue, no flanking residues) with the modifications in square brackets.
//   - A PTM is written as its Unimod accession, [UNIMOD:35], when it has a
//     Unimod id, otherwise as its mass, [+57.0215] (TopPIC's PTM names are
//     not all Unimod names, so they are not written as names).
//   - The protein N-terminal acetylations (M_ACETYLATION, NME_ACETYLATION)
//     and the truncation N-terminal modification are N-terminal
//     modifications: [UNIMOD:1]-SEQ.
//   - Fixed PTMs, the PTMs written in the FASTA sequence and the variable
//     PTMs follow their residue when they are localized to one: C[UNIMOD:4].
//   - Unexpected mass shifts are masses with 4 decimals (sign included) on
//     their residue, or on their range of residues in parentheses,
//     (PEP)[+12.3456]; variable PTMs not localized to one residue likewise,
//     with all the PTMs of the range: (MAM)[UNIMOD:35][UNIMOD:35].
//   - A PTM localized by MIScore (-l) is written at its best site with the
//     site scores of the positive sites (Level 2): S[UNIMOD:21#s1(0.900)]
//     T[#s1(0.100)].
//   - A range that would overlap a range already written (ProForma does
//     not allow it) is written as a set of positions instead, the first
//     residue preferred: P[+12.3456#g1]E[#g1]P[#g1].
//   - A modification without a residue (an empty range) is written as of
//     unknown position: [+12.3456]?SEQ.
// Everything else is Level 1 ProForma.
std::string getProFormaStr(const ProteoformPtr& form_ptr);

}  // namespace proteoform_util

}  // namespace toppic

#endif
