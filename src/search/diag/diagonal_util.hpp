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

#ifndef TOPPIC_SEARCH_DIAG_DIAGONAL_UTIL_HPP_
#define TOPPIC_SEARCH_DIAG_DIAGONAL_UTIL_HPP_

#include "ms/spec/extend_ms.hpp"
#include "prsm/peak_ion_pair.hpp"
#include "search/diag/diag_header.hpp"

namespace toppic {

namespace diagonal_util {

// The precursor refinement of TopMG's graph alignment: the
// precursor-dependent ions paired with peaks within the fragment tolerance
// plus the precursor tolerance T, their errors binned over +-2T, the
// precursor moved to the median bin (+2T with fewer than two pairs).
double refinePrecursorAndHeaderShift(const ProteoformPtr& proteo_ptr,
                                     const ExtendMsPtrVec& ms_three_ptr_vec,
                                     DiagHeaderPtrVec& header_ptrs, double ppo,
                                     double min_mass,
                                     double refine_prec_step_width);

// The precursor refinement of toppic's shift searches (PsAlign,
// PsAlignMultiType). The alignment matched the precursor-dependent ions,
// the C-terminal ions of every diagonal but the last and the N-terminal ions
// of the last, with the fragment tolerance plus the precursor's error
// (N-strict-C-relax and N-relax-C-strict PRM tolerances); the PrSM matches
// them with the fragment tolerance alone. The precursor is moved by the
// delta in [-T, T] (T the precursor tolerance, ppo x the precursor mass)
// that matches the most dependent ions within the fragment tolerance: of
// the best delta interval the one closest to 0, in it the median error of
// its ions; 0 when no ion matches. The headers' shifts are moved as in
// refinePrecursorAndHeaderShift. Returns the refined precursor mass.
double refinePrecursorByMatchedIons(const ProteoformPtr& proteo_ptr,
                                    const ExtendMsPtrVec& ms_three_ptr_vec,
                                    DiagHeaderPtrVec& header_ptrs, double ppo,
                                    double min_mass);

// The precursor refinement of the variable PTM search. Once the alignment
// has chosen the variable PTMs, the proteoform mass is exact: the refined
// precursor mass is the proteoform mass plus ptm_mass, the mass of the
// variable PTMs of the C-terminal diagonal, not a median of fragment errors
// as in refinePrecursorAndHeaderShift. The headers' C-terminal shifts, which
// VarPtmSlowMatch computed from the observed precursor mass, are taken from
// the refined one, so that the C-terminal ions of refineVarPtmHeadersBgnEnd
// carry no precursor error. Returns the refined precursor mass.
double refineVarPtmPrecursorAndHeaderShift(const ProteoformPtr& proteo_ptr,
                                           const DiagHeaderPtrVec& header_ptrs,
                                           double ptm_mass);
// for ps align
DiagHeaderPtrVec refineHeadersBgnEnd(const ProteoformPtr& proteo_ptr,
                                     const ExtendMsPtrVec& ms_three_ptr_vec,
                                     const DiagHeaderPtrVec& heade_ptrs,
                                     double min_mass);
// for ptm align
DiagHeaderPtrVec refineVarPtmHeadersBgnEnd(
    const ProteoformPtr& proteo_ptr, const ExtendMsPtrVec& ms_three_ptr_vec,
    const DiagHeaderPtrVec& header_ptrs, double min_mass);
// for graph alignment
DiagHeaderPtrVec2D refineHeadersBgnEnd(const ProteoformPtr& proteo_ptr,
                                       const ExtendMsPtrVec& ms_three_ptr_vec,
                                       const DiagHeaderPtrVec2D& header_ptrs_2d,
                                       const DiagHeaderPtrVec& header_ptrs_1d,
                                       double min_mass);

int getNewBgn(const PeakIonPairPtrVec& pair_ptrs);

int getNewEnd(const PeakIonPairPtrVec& pair_ptrs);

}  // namespace diagonal_util

} /* namespace toppic */

#endif /* DIAGONAL_HPP_ */
