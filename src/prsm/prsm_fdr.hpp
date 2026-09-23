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

// PrSMs (sorted by E-value) grouped by proteoform cluster, groups ordered by
// their best E-value
PrsmStrPtrVec2D getProteoformGroups(const PrsmStrPtrVec& prsm_ptrs);

// proteoform groups (ordered by best E-value) grouped by protein cluster
PrsmStrPtrVec2D getProteinGroups(const PrsmStrPtrVec2D& proteoforms);

// the proteoform groups of getProteoformGroups (targets and decoys grouped together)
// split by their best PrSM: a proteoform is a target or a decoy according
// to that PrSM
void splitTargetDecoyProteoforms(const PrsmStrPtrVec2D& all_proteoforms,
                                 PrsmStrPtrVec2D& target_proteoforms,
                                 PrsmStrPtrVec2D& decoy_proteoforms);

// spectrum-, proteoform- and protein-level FDRs; the inputs are sorted /
// ordered by E-value
void computePrsmFdr(PrsmStrPtrVec& target_prsms, PrsmStrPtrVec& decoy_prsms);
void computeProteoformFdr(PrsmStrPtrVec2D& target_proteoforms,
                          PrsmStrPtrVec2D& decoy_proteoforms);
void computeProteinFdr(PrsmStrPtrVec2D& target_proteins,
                       PrsmStrPtrVec2D& decoy_proteins);

// Writes every PrSM (targets and decoys) with its spectrum-, proteoform- and
// protein-level FDRs to <base>.<output_file_ext>; the decoys are removed
// later (prsm_decoy_filter) when they are not wanted.
void process(const std::string& spec_file_name,
             const std::string& input_file_ext,
             const std::string& output_file_ext);

}  // namespace prsm_fdr

}  // namespace toppic

#endif
