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
#include "para/peak_tolerance.hpp"
#include "prsm/peak_ion_pair.hpp"
#include "search/diag/diag_header.hpp"

namespace toppic {

namespace diagonal_util {

// The precursor refinement of the shift searches: toppic's one-PTM and
// multiple-PTM searches (PsAlign, PsAlignMultiType) and TopMG's graph
// alignment (GraphAlign, in topmg and toppic --topmg). The alignment matched
// the precursor-dependent ions, the C-terminal ions of every diagonal but the
// last and the N-terminal ions of the last, with the fragment tolerance plus
// the precursor's error: PsAlign and PsAlignMultiType with the
// N-strict-C-relax and N-relax-C-strict PRM tolerances, the graph alignment
// through the distances of its spectrum graph (PRM peaks, a C-terminal one
// P - the fragment mass) within the alignment error tolerance. The PrSM
// matches them with the fragment tolerance alone. The precursor is moved by
// the delta in [-T, T] (T the precursor tolerance: ppo x the precursor mass,
// or the Dalton tolerance of toppic's --prec-dalton-error) that matches the
// most dependent ions within the fragment tolerance: of the best delta
// intervals the one closest to 0, in it the median error of its ions; 0 when
// no ion matches. The last diagonal's N-terminal shift and the other
// diagonals' C-terminal shifts are moved by the delta. A single diagonal
// (only TopMG's graph alignment passes one: an alignment without shifts) is
// anchored at its N-terminus: its C-terminal ions are the dependent ions and
// its C-terminal shift is moved. Returns the refined precursor mass.
double refinePrecursorByMatchedIons(const ProteoformPtr& proteo_ptr,
                                    const ExtendMsPtrVec& ms_three_ptr_vec,
                                    DiagHeaderPtrVec& header_ptrs,
                                    const PeakTolerancePtr& prec_tole_ptr,
                                    double min_mass);

// The precursor refinement of the variable PTM search. Once the alignment
// has chosen the variable PTMs, the proteoform mass is exact: the refined
// precursor mass is the proteoform mass plus ptm_mass, the mass of the
// variable PTMs of the C-terminal diagonal, not a fit to fragment errors as
// in refinePrecursorByMatchedIons. The headers' C-terminal shifts, which
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
