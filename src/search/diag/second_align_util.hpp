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

#ifndef TOPPIC_SEARCH_DIAG_SECOND_ALIGN_UTIL_HPP_
#define TOPPIC_SEARCH_DIAG_SECOND_ALIGN_UTIL_HPP_

#include <string>
#include <vector>

#include "common/base/proteoform_type.hpp"
#include "ms/spec/deconv_ms.hpp"
#include "ms/spec/extend_ms.hpp"
#include "ms/spec/prm_peak.hpp"
#include "para/prsm_para.hpp"
#include "prsm/prsm.hpp"
#include "search/diag/diag_header.hpp"
#include "search/diag/diagonal.hpp"
#include "seq/proteoform.hpp"

namespace toppic {

// The second alignment of toppic's shift searches (hidden --second-align,
// issue #614). The first alignment (PsAlign, PsAlignMultiType) matches the
// precursor-dependent ions with the fragment tolerance plus the precursor's
// (the N-strict-C-relax and N-relax-C-strict PRM tolerances); the PrSM
// matches them with the fragment tolerance at the refined precursor. The
// second alignment repeats the dynamic program at the refined precursor of
// the first alignment's result, on PRM peaks rebuilt at that precursor,
// with the same diagonals moved to it, and (mode strict) with the fragment
// tolerance on every peak, so that the alignment and the annotation agree.
// Its result goes through the same precursor refinement, range refinement
// and PrSM construction as the first one. With no second alignment of a
// positive score, the first alignment's PrSM is kept.
namespace second_align_util {

// main's result of an alignment (PsAlign::geneResult before this branch):
// the precursor refined by the matched ions, the extended spectra at the
// refined precursor, the headers' ranges refined, the PrSM. refine_prec_mass
// and refine_ms_ptr_vec receive the refined precursor and the extended
// spectra at it; nullptr when the range refinement leaves no header.
// With clamp_prec_mass > 0 (the observed precursor mass, for the second
// pass with PrsmPara::useSecondAlignClamp) the refined precursor is kept
// within the precursor tolerance of that mass.
PrsmPtr genePrsm(const ProteoformPtr& proteo_ptr,
                 const DiagHeaderPtrVec& header_ptrs,
                 DeconvMsPtrVec& deconv_ms_ptr_vec,
                 const ExtendMsPtrVec& ms_three_ptr_vec,
                 const PrsmParaPtr& prsm_para_ptr, double& refine_prec_mass,
                 ExtendMsPtrVec& refine_ms_ptr_vec,
                 double clamp_prec_mass = 0);

// The PRM peaks of the spectra at the refined precursor mass (label-free),
// with the zero and precursor peaks, as prm_ms_util::getPrmPeakPtrs builds
// them for the first alignment. With strict, every fragment peak's
// N-strict-C-relax and N-relax-C-strict tolerances are set to its strict
// tolerance (ppo x the fragment mass), the tolerance of the PrSM's matching;
// without it they stay the first alignment's (the fragment tolerance plus
// the precursor's error).
PrmPeakPtrVec genePrmPeaks(const DeconvMsPtrVec& deconv_ms_ptr_vec,
                           const SpParaPtr& sp_para_ptr,
                           double refine_prec_mass, bool strict);

// The first alignment's diagonal headers moved from the observed precursor
// to the refined one (delta = refined - observed): a header anchored at the
// N-terminus (n_strict) keeps its N-terminal shift, one anchored at the
// C-terminus keeps its C-terminal shift, so its N-terminal shift moves by
// delta; every C-terminal shift and the truncation positions are recomputed
// from the refined precursor (DiagHeader::initHeader).
DiagHeaderPtrVec moveHeaders(const DiagonalPtrVec& diagonal_ptrs, double delta,
                             double refine_prec_mass,
                             const ProteoformPtr& proteo_ptr);

// The diagonals of the second alignment (the moved headers' pairs on the
// PRM peaks at the refined precursor, headers without a pair dropped) and
// the peak masses of the dynamic program.
DiagonalPtrVec geneDiagonals(const DiagHeaderPtrVec& header_ptrs,
                             const PrmPeakPtrVec& prm_peaks,
                             int group_spec_num,
                             const ProteoformPtr& proteo_ptr,
                             std::vector<double>& ms_masses);

// Research log of the two passes (environment variable
// TOPPIC_SECOND_ALIGN_LOG = file name): one tab-separated line per
// geneResult call with a first alignment: stage, spectrum id, type, shift
// number, protein, observed precursor, first score, first refined
// precursor, first PrSM (start, end, matched theoretical fragments,
// shifts), second score (- without a second alignment), second refined
// precursor, second PrSM, and whether the first PrSM was kept.
bool logEnabled();

std::string prsmInfo(const PrsmPtr& prsm_ptr);

void log(const std::string& line);

}  // namespace second_align_util

}  // namespace toppic

#endif
