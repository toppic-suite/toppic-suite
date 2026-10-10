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

#ifndef TOPPIC_PRSM_PRSM_FDR_HPP_
#define TOPPIC_PRSM_PRSM_FDR_HPP_

#include <string>

#include "prsm/prsm_str.hpp"

namespace toppic {

namespace prsm_fdr {

// The helpers below are shared with prsm_fdr_groups.

inline bool isDecoy(const PrsmStrPtr& prsm_ptr) {
  return prsm_ptr->getSeqName().find("DECOY_") == 0;
}

// The score the FDRs rank PrSMs by, smaller is better: the E-value, or with
// --rescore the rescoring step's score negated (prsm/prsm_rescore).
using FdrKey = double (*)(const PrsmStrPtr&);

inline double eValueKey(const PrsmStrPtr& prsm_ptr) {
  return prsm_ptr->getEValue();
}

inline double rescoreKey(const PrsmStrPtr& prsm_ptr) {
  return -prsm_ptr->getRescore();
}

// PrSMs (sorted by E-value) grouped by proteoform cluster, groups ordered by
// their best E-value
PrsmStrPtrVec2D getProteoformGroups(const PrsmStrPtrVec& prsm_ptrs);

// proteoform groups (ordered by best E-value) grouped by protein cluster
PrsmStrPtrVec2D getProteinGroups(const PrsmStrPtrVec2D& proteoforms);

// the proteoform groups of getProteoformGroups (targets and decoys grouped
// together) split by their best PrSM: a proteoform is a target or a decoy
// according to that PrSM
void splitTargetDecoyProteoforms(const PrsmStrPtrVec2D& all_proteoforms,
                                 PrsmStrPtrVec2D& target_proteoforms,
                                 PrsmStrPtrVec2D& decoy_proteoforms);

// spectrum-, proteoform- and protein-level FDRs; the inputs are sorted /
// ordered by the key (the E-value by default)
void computePrsmFdr(PrsmStrPtrVec& target_prsms, PrsmStrPtrVec& decoy_prsms,
                    FdrKey key = eValueKey);
void computeProteoformFdr(PrsmStrPtrVec2D& target_proteoforms,
                          PrsmStrPtrVec2D& decoy_proteoforms,
                          FdrKey key = eValueKey);
void computeProteinFdr(PrsmStrPtrVec2D& target_proteins,
                       PrsmStrPtrVec2D& decoy_proteins,
                       FdrKey key = eValueKey);

// The FDR cutoffs at which the two rounds of a rescored FDR computation
// count the targets that pass (spectrum, proteoform and protein level).
struct FdrCutoffs {
  double prsm = 0.01;
  double proteoform = 0.01;
  double protein = 0.01;
};

// Writes every PrSM (targets and decoys) with its spectrum-, proteoform- and
// protein-level FDRs to <base>.<output_file_ext>; the decoys are removed
// later (prsm_decoy_filter) when they are not wanted. use_rescore: the
// PrSMs are ranked by the rescoring step's score (default with -d,
// prsm/prsm_rescore) instead of the E-value; every PrSM must have one.
// rescore_fallback (with use_rescore; off with the hidden
// --disable-rescore-fallback): two rounds, the FDRs by the E-value and by
// the rescore, and at each level (PrSM, proteoform, protein) the rescored
// FDRs unless the E-value's give more targets within the level's cutoff
// (toppic_test experiment 153: the rescoring raised the Bruker run's PrSMs
// but lowered its proteoforms and proteins, and collapsed on a data set
// with few decoys); a tie keeps the rescore. Each search (toppic, topmg)
// chooses on its own.
void process(const std::string& spec_file_name,
             const std::string& input_file_ext,
             const std::string& output_file_ext, bool use_rescore = false,
             const FdrCutoffs& cutoffs = FdrCutoffs(),
             bool rescore_fallback = true);

}  // namespace prsm_fdr

}  // namespace toppic

#endif
