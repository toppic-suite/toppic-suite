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

#ifndef TOPPIC_FILTER_MASSMATCH_MASS_MATCH_UTIL_HPP_
#define TOPPIC_FILTER_MASSMATCH_MASS_MATCH_UTIL_HPP_

#include "filter/massmatch/mass_match.hpp"
#include "seq/prot_candidate.hpp"

namespace toppic {

namespace mass_match_util {

ProtCandidatePtrVec findZeroShiftTopProteins(std::vector<short>& scores,
                                             std::vector<short>& rev_scores,
                                             const MassMatchPtr& index_ptr,
                                             const MassMatchPtr& rev_index_ptr,
                                             double prec_minus_water_mass,
                                             double prec_error_tole,
                                             int threshold, int num);

// the simple version is faster, but less accurate
// max_scores[i]: the best row score of proteoform i in the index
void compProteoMaxScores(const std::vector<short>& scores,
                         const MassMatchPtr& index_ptr,
                         std::vector<int>& max_scores);

// Zero-shift candidates: the proteoforms whose best N-terminal row score
// plus best C-terminal row score reaches the threshold (no truncation shift
// is recorded). The maxima come from compProteoMaxScores, so the four
// proteoform types share four scans of the score arrays.
ProtCandidatePtrVec simpleFindZeroShiftTopProteins(
    const std::vector<int>& max_scores, const std::vector<int>& rev_max_scores,
    int threshold, int top_num);

// One-shift candidates: the proteoforms with the best total score of an
// N-terminal row and a C-terminal row whose masses leave a shift in
// [min_shift, max_shift] of the precursor mass. The per-proteoform maxima
// come from compProteoMaxScores, so the four proteoform types share four
// scans of the score arrays, and they bound the total score of a proteoform,
// which keeps the truncation row scan to the proteoforms that can be
// reported.
ProtCandidatePtrVec findOneShiftTopProteins(
    std::vector<short>& scores, std::vector<short>& rev_scores,
    const std::vector<int>& max_scores, const std::vector<int>& rev_max_scores,
    const MassMatchPtr& index_ptr, const MassMatchPtr& rev_index_ptr,
    double prec_minus_water_mass, double prec_error_tole, double min_shift,
    double max_shift, int term_shift_num, int threshold, int top_num);

ProtCandidatePtrVec findVarPtmTopProteins(
    std::vector<short>& scores, std::vector<short>& rev_scores,
    const MassMatchPtr& index_ptr, const MassMatchPtr& rev_index_ptr,
    double prec_minus_water_mass, double prec_error_tole,
    std::vector<double>& ptm_shifts, int threshold, int top_num);

// prec_match_rows[i]: the rows of proteoform i that received the precursor
// bonus, in ascending row order, collected from the rows the scoring noted
// instead of a pass over every row of the index. A row without the bonus
// scores at most the number of masses of the spectrum, far below the bonus,
// and the threshold of the variable-PTM filter is twice the bonus, so no
// other row can belong to a candidate.
void compPrecMatchRows(const MassMatchPtr& index_ptr,
                       const std::vector<int>& prec_rows,
                       std::vector<std::vector<int>>& prec_match_rows);

// The same candidates from those row lists, which the four proteoform types
// of a spectrum share.
ProtCandidatePtrVec findVarPtmTopProteins(
    const std::vector<short>& scores, const std::vector<short>& rev_scores,
    const std::vector<std::vector<int>>& prec_match_rows,
    const std::vector<std::vector<int>>& rev_prec_match_rows,
    const MassMatchPtr& index_ptr, const MassMatchPtr& rev_index_ptr,
    double prec_minus_water_mass, double prec_error_tole,
    std::vector<double>& ptm_shifts, int threshold, int top_num);

ProtCandidatePtrVec findDiagTopProteins(
    const std::vector<short>& scores, const std::vector<int>& proteo_row_begins,
    const std::vector<int>& proteo_row_ends, int threshold, int top_num);

}  // namespace mass_match_util

}  // namespace toppic

#endif
