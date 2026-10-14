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

#ifndef TOPPIC_STAT_LOCAL_ALIGN_LOCAL_ALIGN_PROCESSOR_HPP_
#define TOPPIC_STAT_LOCAL_ALIGN_LOCAL_ALIGN_PROCESSOR_HPP_

#include <fstream>
#include <string>

#include "prsm/prsm.hpp"
#include "stat/local/local_mng.hpp"
#include "stat/local_align/local_align_util.hpp"

namespace toppic {

// The alternative to LocalProcessor (toppic --local-align with -B): the
// unexpected shift of a PrSM is explained by the PTMs of the -B file, one
// PTM or an unordered pair, within the precursor tolerance; it is replaced
// by the PTM(s), as variable PTMs with the MIScore of LocalProcessor, only
// when exactly one explanation exists, and without LocalProcessor's filters
// on the matched fragments. A PrSM with two shifts is split at the middle
// matched cleavage site between them and each piece is explained on its
// own.
class LocalAlignProcessor {
 public:
  explicit LocalAlignProcessor(const LocalMngPtr& mng_ptr);

  void process();

 private:
  void init();

  // the switches of TOPPIC_LOCAL_ALIGN_PARA ("name=value,...")
  void readPara();

  // The result of a piece: the new VARIABLE shifts, in the coordinates of
  // cand_form when the proteoform's start or end was moved to explain the
  // shift (nullptr: the original proteoform's)
  struct PieceResult {
    bool ok = false;
    MassShiftPtrVec shifts;
    ProteoformPtr cand_form;
    std::string decision;
    // the other explanations of the shift on the same proteoform as the
    // chosen one (cand_form; the proteoform as it is when none was
    // chosen), each localized: Proteoform::alternative_mass_shift_list_
    std::vector<MassShiftPtrVec> alternatives;
  };

  // A candidate explanation of a shift: the proteoform it is on (nullptr:
  // the proteoform as it is), the piece on that proteoform, the
  // explanation, whether the proteoform carries the N-terminal
  // acetylation, and the name for the log
  struct Cand {
    ProteoformPtr form;
    local_align::Piece piece;
    local_align::Explanation expl;
    bool acetyl_form = false;
    std::string name;
  };

  // The tie-break rules (local_align::breakTies) and, for the explanations
  // left that differ only by Acetyl vs Trimethyl on the same proteoform,
  // the site scores (both localized, the better match score wins, then
  // the higher confidence, then Acetyl), then the one-explanation rule.
  // Returns the index of the chosen candidate with its localized shifts
  // and match score, or -1 (ambiguous); rule names the rules that decided.
  int chooseOne(const PrsmPtr& prsm_ptr, const std::vector<Cand>& cands,
                MassShiftPtrVec& shifts, int& best_score, std::string& rule);

  // chooseOne, plus the other candidates on the chosen candidate's
  // proteoform (the proteoform as it is when none is chosen) localized as
  // the alternatives
  int choose(const PrsmPtr& prsm_ptr, const std::vector<Cand>& cands,
             MassShiftPtrVec& shifts, int& best_score, std::string& rule,
             std::vector<MassShiftPtrVec>& alternatives);

  PrsmPtr processOnePrsm(const PrsmPtr& prsm_ptr);

  // The alternative mass shift lists of the result proteoform out_form
  // (Proteoform::alternative_mass_shift_list_): for every alternative of a
  // shift, the proteoform's shifts other than the unexpected ones and the
  // variable PTMs the localization added (new_shifts), the alternative's
  // variable PTMs and the other unexpected shift(s) as they are. ori_start
  // and unexp_shifts are the original proteoform's.
  void setAlternatives(
      const ProteoformPtr& out_form, int ori_start,
      const MassShiftPtrVec& unexp_shifts, const MassShiftPtrVec& new_shifts,
      const std::vector<std::pair<MassShiftPtr, PieceResult>>& alt_results);

  // The shift of the piece explained and localized on the proteoform as it
  // is and, with n_forms (the piece holds the N-terminus and the shift has
  // few matched sites on its N-terminal side), on its other N-terminal
  // forms of the same start (the acetylation dropped or added); the
  // tie-break rules over all candidates, then the one-explanation rule.
  // kept_shifts: the proteoform's shifts other than the one explained.
  PieceResult processPiece(const PrsmPtr& prsm_ptr,
                           const local_align::Piece& piece,
                           const std::string& piece_name, int shift_idx,
                           const MassShiftPtrVec& kept_shifts, bool n_forms);

  // The shift explained on the candidate proteoforms whose start (n_side)
  // or end is moved by up to term_range residues: the tie-break rules,
  // then the one-explanation rule over all candidates. piece_bgn /
  // piece_end are the piece's bounds in the original proteoform,
  // kept_shifts its shifts other than the one explained (and the protein
  // N-terminal ones).
  PieceResult processTerminal(const PrsmPtr& prsm_ptr,
                              const MassShiftPtr& shift_ptr, bool n_side,
                              int piece_bgn, int piece_end,
                              const MassShiftPtrVec& kept_shifts,
                              const std::string& piece_name, int shift_idx);

  // rule: the tie-break rules that chose the explanation (breakTies);
  // all_names: the explanations before the tie-break when it applied
  void log(const PrsmPtr& prsm_ptr, int shift_idx, double shift_mass,
           const std::string& piece_name, const std::string& expl_names,
           const std::string& decision, const std::string& sites,
           const std::string& rule = "", const std::string& all_names = "");

  LocalMngPtr mng_ptr_;

  local_align::Para para_;

  std::ofstream log_file_;
};

using LocalAlignProcessorPtr = std::shared_ptr<LocalAlignProcessor>;

}  // namespace toppic

#endif
