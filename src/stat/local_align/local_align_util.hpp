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

#ifndef TOPPIC_STAT_LOCAL_ALIGN_LOCAL_ALIGN_UTIL_HPP_
#define TOPPIC_STAT_LOCAL_ALIGN_LOCAL_ALIGN_UTIL_HPP_

#include <string>
#include <vector>

#include "common/base/mod.hpp"
#include "common/base/ptm.hpp"
#include "ms/spec/extend_ms.hpp"
#include "seq/mass_shift.hpp"
#include "seq/proteoform.hpp"
#include "stat/local/local_mng.hpp"

namespace toppic {

namespace local_align {

// The PTM list of -B and the switches of the localization
struct Para {
  // the PTMs of the -B file, and its modifications by position
  PtmPtrVec ptm_ptr_vec;
  ModPtrVec mod_list_N;
  ModPtrVec mod_list_C;
  ModPtrVec mod_list_any;
  // a PTM mass may match the shift at one isotope mass away (off: a 1 Da
  // precursor error is not turned into a PTM; toppic_test experiment 140)
  bool iso = false;
  // the explanations have at most this many PTMs (1 or 2)
  int max_ptm = 2;
  // two-shift PrSMs are split and explained piece by piece
  bool split = true;
  // a shift without an explanation that has at most term_max_matched
  // matched cleavage sites on its N-terminal (C-terminal) side is tried on
  // the proteoforms whose start (end) is moved by up to term_range residues,
  // and on the N-terminal side also on the other N-terminal forms of the
  // unmoved start (the N-terminal acetylation dropped or added)
  bool term = true;
  int term_range = 5;
  int term_max_matched = 3;
  // several explanations are reduced by the tie-break rules of breakTies
  // before the one-explanation rule (0 = the pure rule: several
  // explanations always stop the localization)
  bool tie = true;
};

// A shift of a piece: its cleavage-site range [left, right) in the piece's
// positions and its mass
struct ShiftRange {
  int left = 0;
  int right = 0;
  double mass = 0;
};

// A proteoform, or the part of one between two unexpected shifts, with one
// unexpected shift to explain. The positions are those of the proteoform's
// residues; the theoretical masses of the piece's cleavage sites carry the
// piece's known shifts and the masses outside the piece, so that they are
// matched against the real spectra.
struct Piece {
  // the first residue of the piece in the proteoform and the one after the
  // last; the proteoform's length
  int bgn = 0;
  int end = 0;
  int form_len = 0;
  ResiduePtrVec residues;
  // the mass of the unexpected shift to explain
  double shift_mass = 0;
  // the variable and the other unexpected shifts inside the piece, and the
  // fixed ones, with the piece's positions
  std::vector<ShiftRange> known_shifts;
  std::vector<ShiftRange> fixed_shifts;
  // the residue masses and shifts before (after) the piece: added to the
  // prefix (suffix) fragment masses
  double n_offset = 0;
  double c_offset = 0;
  // cleavage site j (0..len): the residue masses and known shifts before it
  // (prefix) and from it on (suffix); a site inside a known shift's range
  // gives no fragment
  std::vector<double> prefix_mass;
  std::vector<double> suffix_mass;
  std::vector<bool> valid;

  int getLen() const { return end - bgn; }
};

// An explanation of a shift: one PTM, or two (second != nullptr)
struct Explanation {
  PtmPtr first;
  PtmPtr second;

  std::string getName() const;
  int getPtmNum() const { return second == nullptr ? 1 : 2; }
  // a PTM of the explanation has this abbreviation (the -B file's name)
  bool hasPtm(const std::string& abbr) const;
};

// The tie-break rules for several explanations of one shift, applied in
// this order while more than one is left: (1) fewer PTMs (Dimethyl over
// Methyl+Methyl, Trimethyl over Methyl+Dimethyl, a PTM on the N-terminal
// acetylation form over Acetyl + the PTM on the plain form); (2) the
// N-terminal acetylation form over an Acetyl of -B on the plain form
// (acetyl_forms[i]: the explanation's candidate proteoform carries the
// acetylation). Acetyl vs Trimethyl (0.036 Da apart) is not decided here:
// LocalAlignProcessor::choose localizes both and takes the better match
// score. Returns the indices of the explanations kept; rule names the
// rules that reduced the set ("fewer_ptms", "acetyl_form", joined by '+';
// empty when none applied).
std::vector<size_t> breakTies(const std::vector<Explanation>& expls,
                              const std::vector<bool>& acetyl_forms,
                              std::string& rule);

// the shifts of the proteoform that count as known in the pieces
bool isKnownShift(const MassShiftPtr& shift_ptr);

// The piece of the residues [bgn, end) of the proteoform that explains
// shift_ptr; the shifts outside the piece give its offsets.
Piece makePiece(const ProteoformPtr& form_ptr, int bgn, int end,
                const MassShiftPtr& shift_ptr);

// The cleavage sites of the proteoform (with all its shifts) that a
// fragment of the spectra matches, as a 0/1 vector over the sites 0..len
std::vector<int> compMatchedSites(const ProteoformPtr& form_ptr,
                                  const ExtendMsPtrVec& extend_ms_ptr_vec,
                                  const LocalMngPtr& mng_ptr);

// the number of matched sites in [from, to]
int countMatchedSites(const std::vector<int>& matched_sites, int from, int to);

// the shift moved by -delta positions, with its local annotation
MassShiftPtr shiftPositions(const MassShiftPtr& shift_ptr, int delta);

// The break site between the two unexpected shifts: the middle matched site
// of [first shift's right position, second shift's left position], or the
// middle of the range when none is matched. matched_num gets the number of
// matched sites in the range.
int findBreakSite(const MassShiftPtr& shift_1, const MassShiftPtr& shift_2,
                  const std::vector<int>& matched_sites, int& matched_num);

bool modifiable(const Piece& piece, int i, const PtmPtr& ptm_ptr,
                const Para& para);

// the explanations of the piece's shift: the single PTMs and unordered PTM
// pairs of the list whose mass matches the shift and that have sites
std::vector<Explanation> findExplanations(const Piece& piece, double err_tole,
                                          const Para& para);

// The shift of the piece replaced by the explanation: one or two VARIABLE
// mass shifts with local annotations, in the proteoform's positions.
// best_score gets the matched fragment score of the localization.
MassShiftPtrVec localize(const Piece& piece, const Explanation& expl,
                         const ExtendMsPtrVec& extend_ms_ptr_vec,
                         const Para& para, const LocalMngPtr& mng_ptr,
                         int& best_score);

}  // namespace local_align

}  // namespace toppic

#endif
