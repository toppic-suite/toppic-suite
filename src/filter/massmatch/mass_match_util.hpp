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

#include <cstdint>
#include <vector>

#include "common/base/prot_mod.hpp"
#include "filter/massmatch/mass_match.hpp"
#include "seq/prot_candidate.hpp"

namespace toppic {

namespace mass_match_util {

// Whether an index read from a file was built at a filter's column limits
// (the largest proteoform mass and the scale of the integer masses).
bool hasFilterLimits(const MassMatchPtr& index_ptr, int max_proteoform_mass,
                     int filter_scale);

// Prints, once per run, that the index file does not have the filter's
// limits and that the filter builds its indexes in memory instead.
void warnIndexLimits(const std::string& file_name,
                     const MassMatchPtr& index_ptr, int max_proteoform_mass,
                     int filter_scale);

// Prints, once per run, that the index directory holds the files of an
// earlier version, which wrote the diagonal index the one-shift and
// multiple-shift filters share as one_ptm_diag_index and multi_ptm_index
// instead of ptm_diag_index (suffix: the parameter and block part of the
// file names), so that those filters build their indexes in memory.
// Prints nothing when neither old file exists.
void warnOldIndexLayout(const std::string& index_dir,
                        const std::string& suffix);

// The sequence ranges of the index rows (issue #348). A pair of a row of an
// N-terminal side index and a row of a C-terminal side index is a range of
// the protein, [start, end], 0-based positions of its first and last
// residue.
//
// The first residue of each row of an N-terminal side index: a row of the
// N-terminal form index (is_term_index) starts at the truncation length of
// its protein modification, found by its shift (0 for NONE and
// M_ACETYLATION, 1 for NME and NME_ACETYLATION); the n-th row of a
// proteoform in the diagonal index removes its first n residues and starts
// at n.
std::vector<int> compRowStartPositions(const MassMatchPtr& index_ptr,
                                       bool is_term_index,
                                       const ProtModPtrVec& prot_mod_ptrs);

// The last residue of each row of a C-terminal side index: the one row of a
// proteoform in the reverse term index (is_term_index) ends at the last
// residue, length - 1; the n-th row in the reverse diagonal index removes
// the last n residues and ends at length - 1 - n.
std::vector<int> compRowEndPositions(const MassMatchPtr& index_ptr,
                                     bool is_term_index,
                                     const std::vector<int>& proteo_lens);

// The ranges of the rows of a filter's two indexes for one proteoform type,
// and the ranges the type allows: with proper_start the range may not start
// at the first residue (a zero-shift SUFFIX or INTERNAL proteoform), with
// proper_end it may not end at the last one (PREFIX, INTERNAL).
struct RangeRows {
  const std::vector<int>* row_starts_ = nullptr;
  const std::vector<int>* rev_row_ends_ = nullptr;
  const std::vector<int>* proteo_lens_ = nullptr;
  bool proper_start_ = false;
  bool proper_end_ = false;
};

// The row ranges of a filter's four indexes: the N-terminal form index, the
// diagonal index, the reverse term index and the reverse diagonal index,
// with the lengths of the proteins, computed once per database block.
struct IndexRowRanges {
  std::vector<int> proteo_lens_;
  std::vector<int> term_starts_;
  std::vector<int> diag_starts_;
  std::vector<int> rev_term_ends_;
  std::vector<int> rev_diag_ends_;

  void init(const std::vector<int>& proteo_lens, const MassMatchPtr& term_ptr,
            const MassMatchPtr& diag_ptr, const MassMatchPtr& rev_term_ptr,
            const MassMatchPtr& rev_diag_ptr,
            const ProtModPtrVec& prot_mod_ptrs);

  // The rows of a proteoform type: the N-terminal form index (term_start)
  // or the diagonal index, the reverse term index (term_end) or the reverse
  // diagonal index.
  RangeRows getRangeRows(bool term_start, bool term_end, bool proper_start,
                         bool proper_end) const;
};

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

// One-shift candidates: the proteoforms with the best total score of an
// N-terminal row and a C-terminal row whose masses leave a shift in
// [min_shift, max_shift] of the precursor mass. The per-proteoform maxima
// come from compProteoMaxScores, so the four proteoform types share four
// scans of the score arrays, and they bound the total score of a proteoform,
// which keeps the truncation row scan to the proteoforms that can be
// reported.
// For each proteoform of an index, whether its truncation shifts do not
// increase from one row to the next (1) or do somewhere (0).
std::vector<char> compNonincreasingShifts(const MassMatchPtr& index_ptr);

// nonincr_shifts, when given, is compNonincreasingShifts of index_ptr: the
// truncation window of a proteoform whose shifts do not increase is then
// moved only for the rows that are scanned, which gives the same window.
// range_rows, when given, sets each candidate's range to that of its best
// pair.
ProtCandidatePtrVec findOneShiftTopProteins(
    std::vector<short>& scores, std::vector<short>& rev_scores,
    const std::vector<int>& max_scores, const std::vector<int>& rev_max_scores,
    const MassMatchPtr& index_ptr, const MassMatchPtr& rev_index_ptr,
    double prec_minus_water_mass, double prec_error_tole, double min_shift,
    double max_shift, int term_shift_num, int threshold, int top_num,
    const std::vector<char>* nonincr_shifts = nullptr,
    const RangeRows* range_rows = nullptr);

// Candidates as sequence ranges (issue #348): a candidate is a protein and a
// range, so a protein can be reported once for each of its ranges, and the
// pairs of rows that give the same range (two N-terminal forms with the same
// start) give one candidate, the better. The best top_num candidates are
// reported, by score descending, then protein id, start and end ascending.
//
// One-shift: the pairs whose range mass is in the shift window; a
// candidate's shift lists hold the N- and the C-terminal shift of its own
// pair.
ProtCandidatePtrVec findOneShiftTopRanges(
    std::vector<short>& scores, std::vector<short>& rev_scores,
    const std::vector<int>& max_scores, const std::vector<int>& rev_max_scores,
    const MassMatchPtr& index_ptr, const MassMatchPtr& rev_index_ptr,
    const RangeRows& range_rows, double prec_minus_water_mass,
    double prec_error_tole, double min_shift, double max_shift,
    int term_shift_num, int threshold, int top_num,
    const std::vector<char>* nonincr_shifts = nullptr);

ProtCandidatePtrVec findVarPtmTopProteins(
    std::vector<short>& scores, std::vector<short>& rev_scores,
    const MassMatchPtr& index_ptr, const MassMatchPtr& rev_index_ptr,
    double prec_minus_water_mass, double prec_error_tole,
    std::vector<double>& ptm_shifts, int threshold, int top_num);

// The rows of an index that received the precursor bonus, grouped by
// proteoform: rows_ in ascending order, and the rows of proteoform i are
// rows_[begins_[i]] to rows_[begins_[i + 1] - 1] (the rows of a proteoform
// are one range of the index, so ascending order groups them). words_ is
// scratch space, kept with the lists so that it is allocated once per block.
struct PrecMatchRows {
  std::vector<int> rows_;
  std::vector<int> begins_;
  std::vector<uint64_t> words_;
};

// The rows of each proteoform that received the precursor bonus, collected
// from the rows the scoring noted instead of a pass over every row of the
// index. A row without the bonus scores at most the number of masses of the
// spectrum, far below the bonus, and the threshold of the variable-PTM filter
// is twice the bonus, so no other row can belong to a candidate.
void compPrecMatchRows(const MassMatchPtr& index_ptr,
                       const std::vector<int>& prec_rows,
                       PrecMatchRows& prec_match_rows);

// The same candidates from those row lists, which the four proteoform types
// of a spectrum share. range_rows, when given, sets each candidate's range to
// that of its best pair. start_row_mask, when given, has one entry per row of
// the first index, and only the rows whose entry is nonzero start a pair (the
// rows whose first residue can carry a proteoform N-terminal PTM).
ProtCandidatePtrVec findVarPtmTopProteins(
    const std::vector<short>& scores, const std::vector<short>& rev_scores,
    const PrecMatchRows& prec_match_rows,
    const PrecMatchRows& rev_prec_match_rows,
    const MassMatchPtr& index_ptr, const MassMatchPtr& rev_index_ptr,
    double prec_minus_water_mass, double prec_error_tole,
    std::vector<double>& ptm_shifts, int threshold, int top_num,
    const RangeRows* range_rows = nullptr,
    const std::vector<char>* start_row_mask = nullptr);

// Variable-PTM candidates as sequence ranges (see findOneShiftTopRanges):
// the pairs whose range mass plus one of the PTM shifts is within the
// precursor error tolerance. start_row_mask as in findVarPtmTopProteins.
ProtCandidatePtrVec findVarPtmTopRanges(
    const std::vector<short>& scores, const std::vector<short>& rev_scores,
    const PrecMatchRows& prec_match_rows,
    const PrecMatchRows& rev_prec_match_rows,
    const MassMatchPtr& index_ptr, const MassMatchPtr& rev_index_ptr,
    const RangeRows& range_rows, double prec_minus_water_mass,
    double prec_error_tole, std::vector<double>& ptm_shifts, int threshold,
    int top_num, const std::vector<char>* start_row_mask = nullptr);

// Zero-shift candidates as sequence ranges (see findOneShiftTopRanges): the
// pairs of rows that both received the precursor bonus, whose range mass is
// within the precursor error tolerance and whose range the proteoform type
// allows. With one_per_protein, each protein keeps its best pair, the first
// in the order of the rows, as findZeroShiftTopProteins does.
ProtCandidatePtrVec findZeroShiftTopRanges(
    const std::vector<short>& scores, const std::vector<short>& rev_scores,
    const PrecMatchRows& prec_match_rows,
    const PrecMatchRows& rev_prec_match_rows,
    const MassMatchPtr& index_ptr, const MassMatchPtr& rev_index_ptr,
    const RangeRows& range_rows, double prec_minus_water_mass,
    double prec_error_tole, int threshold, int top_num, bool one_per_protein);

ProtCandidatePtrVec findDiagTopProteins(
    const std::vector<short>& scores, const std::vector<int>& proteo_row_begins,
    const std::vector<int>& proteo_row_ends, int threshold, int top_num);

}  // namespace mass_match_util

}  // namespace toppic

#endif
