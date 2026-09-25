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

#ifndef TOPPIC_SEARCH_PTMSEARCH_PS_ALIGN_MULTI_TYPE_HPP_
#define TOPPIC_SEARCH_PTMSEARCH_PS_ALIGN_MULTI_TYPE_HPP_

#include <array>
#include <vector>

#include "common/base/proteoform_type.hpp"
#include "prsm/prsm.hpp"
#include "search/diag/diag_header.hpp"
#include "search/diag/diagonal.hpp"
#include "search/oneptmsearch/ps_align_para.hpp"
#include "search/ptmsearch/multi_type_dp_pair.hpp"

namespace toppic {

// A faster PsAlign (search/oneptmsearch/ps_align.hpp), for
// PtmSearchSlowFilter / PtmSlowMatch only: the same diagonal-pair dynamic
// program (the same recursion, the same diagonal/shift/truncation
// transitions, the same idxes_/penalties_ shift-candidate precomputation),
// computed for the COMPLETE, PREFIX, SUFFIX and INTERNAL proteoform types
// together in one pass over the diagonal pairs, instead of PsAlign's one
// pass per type.
//
// The two only differ in the two truncation rules of getTruncPre_ (below)
// -- whether the alignment may start/end at a protein terminus or at any
// peptide bond -- which are still applied per type, exactly as
// PsAlign::getTruncPre does (COMPLETE: protein N-term start, protein C-term
// end; PREFIX: protein N-term start, any C-term end; SUFFIX: any N-term
// start, protein C-term end; INTERNAL: any start, any end). The interior
// diagonal and shift transitions do not depend on the type at all and are
// shared, which is almost all of the dynamic program's cost.
//
// This follows the single dynamic program of the original alignment
// algorithm this was ported from,
// tools/proteomics_msalign_plus_0.7-master's
// edu.ucsd.msalign.align.ptm.PSAlign, whose dp() loops the pair position p,
// the shift count s and the type t = 0..3 together in one nest -- while
// keeping PsAlign's position-indexed shift-candidate lists (idxes_/
// penalties_, built once in initDpPair_) in place of PSAlign's O(number of
// earlier pairs) shift scan, and PsAlign's caching of the type-independent
// diagonal predecessor.
class PsAlignMultiType {
 public:
  PsAlignMultiType(const std::vector<double>& ms_masses,
                   const std::vector<double>& seq_masses,
                   const DiagonalPtrVec& diagonal_ptrs,
                   const PsAlignParaPtr& para_ptr);

  // Runs the dynamic program and the backtrace for exactly the given
  // proteoform types (a non-empty subset of COMPLETE/PREFIX/SUFFIX/
  // INTERNAL, any order): only those types' table rows are computed.
  void compute(const ProteoformTypePtrVec& type_ptrs);

  double getAlignScr(const ProteoformTypePtr& type_ptr, int s) {
    return align_scores_[type_ptr->getId()][s];
  }

  DiagHeaderPtrVec getDiagHeaders(const ProteoformTypePtr& type_ptr, int s) {
    return backtrack_diagonal_ptrs_[type_ptr->getId()][s];
  }

  PrsmPtr geneResult(const ProteoformTypePtr& type_ptr, int shift_num,
                     const ProteoformPtr& proteo_ptr,
                     DeconvMsPtrVec& deconv_ms_ptr_vec,
                     ExtendMsPtrVec& ms_three_ptr_vec,
                     const PrsmParaPtr& prsm_para_ptr);

 protected:
  PsAlignParaPtr para_ptr_;

  std::vector<double> ms_masses_;

  std::vector<double> seq_masses_;

  DiagonalPtrVec diagonal_ptrs_;

  std::vector<std::vector<int>> idxes_;

  std::vector<std::vector<bool>> penalties_;

  std::vector<MultiTypeDpPairPtrVec> dp_2d_pair_ptrs_;

  MultiTypeDpPairPtr first_pair_ptr_;

  MultiTypeDpPairPtr last_pair_ptr_;

  MultiTypeDpPairPtrVec segment_bgn_pair_ptrs_;

  MultiTypeDpPairPtrVec segment_end_pair_ptrs_;

  MultiTypeDpPairPtrVec dp_pair_ptrs_;

  // [type_id][s]
  std::vector<DiagHeaderPtrVec2D> backtrack_diagonal_ptrs_;

  std::vector<std::vector<double>> align_scores_;

  void initDpPair();

  void dpPrep(const ProteoformTypePtrVec& type_ptrs);

  void dp(const ProteoformTypePtrVec& type_ptrs);

  void backtrace(const ProteoformTypePtrVec& type_ptrs);

  DiagHeaderPtrVec backtrace(int type_id, int s);

  // Fills trunc_prev/trunc_score, one entry per requested type (indexed by
  // ProteoformType::getId()), with the best segment-end predecessor of the
  // final (virtual) pair for each type and its score at shift count s. Scans
  // segment_end_pair_ptrs_ once for all requested types together, since only
  // the two end-truncation flags read off each candidate (isProtCTermMatch,
  // isPepCTermMatch) are needed to serve every type -- mirrors getShiftPre's
  // fusion of the per-diagonal scan in dp().
  void getTruncPreLast(
      const ProteoformTypePtrVec& type_ptrs, int s,
      std::array<MultiTypeDpPairPtr, kProteoformTypeNum>* trunc_prev,
      std::array<double, kProteoformTypeNum>* trunc_score);
};

using PsAlignMultiTypePtr = std::shared_ptr<PsAlignMultiType>;

} /* namespace toppic */

#endif
