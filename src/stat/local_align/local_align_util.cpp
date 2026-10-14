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

#include "stat/local_align/local_align_util.hpp"

#include <algorithm>
#include <cmath>
#include <functional>

#include "common/base/amino_acid_base.hpp"
#include "common/base/mass_constant.hpp"
#include "common/base/mod_type.hpp"
#include "common/base/residue_base.hpp"
#include "common/util/logger.hpp"
#include "ms/factory/extend_ms_util.hpp"
#include "seq/alter.hpp"
#include "seq/alter_type.hpp"
#include "seq/local_anno.hpp"
#include "stat/local/local_util.hpp"

namespace toppic {

namespace local_align {

namespace {

// the score of a layer of the two-PTM DP before any site is modifiable
const int NEG_SCORE = -100000000;

// 1 when a mass of the sorted spectrum masses is within the tolerance of
// the theoretical mass
int matchMass(double theo_mass, const std::vector<double>& exp_masses,
              const PeakTolerancePtr& tole_ptr) {
  double err = tole_ptr->compStrictErrorTole(theo_mass);
  auto it = std::lower_bound(exp_masses.begin(), exp_masses.end(),
                             theo_mass - err);
  if (it != exp_masses.end() && *it <= theo_mass + err) {
    return 1;
  }
  return 0;
}

// The matches of the cleavage sites 1..len-1 of the piece over the spectra:
// row_n[j] counts the spectra with a prefix fragment of site j (the piece's
// prefix mass plus the N-terminal offset and ion shift plus n_add), row_c[j]
// those with a suffix fragment (the suffix mass plus the C-terminal offset,
// ion shift and water plus c_add).
void compMatchRows(const Piece& piece, const ExtendMsPtrVec& extend_ms_ptr_vec,
                   double n_add, double c_add, const LocalMngPtr& mng_ptr,
                   std::vector<int>& row_n, std::vector<int>& row_c) {
  int len = piece.getLen();
  row_n.assign(len + 1, 0);
  row_c.assign(len + 1, 0);
  PeakTolerancePtr tole_ptr = mng_ptr->frag_tole_ptr_;
  double water = mass_constant::getWaterMass();
  for (size_t s = 0; s < extend_ms_ptr_vec.size(); s++) {
    std::vector<double> ms_masses =
        extend_ms_util::getExtendMassVec(extend_ms_ptr_vec[s]);
    ActivationPtr act_ptr =
        extend_ms_ptr_vec[s]->getMsHeaderPtr()->getActivationPtr();
    double n_shift = act_ptr->getN_BYShift() + piece.n_offset + n_add;
    double c_shift = act_ptr->getC_BYShift() + water + piece.c_offset + c_add;
    // the termini (sites 0 and len) are the precursor, not fragments
    for (int j = 1; j < len; j++) {
      if (!piece.valid[j]) continue;
      row_n[j] += matchMass(piece.prefix_mass[j] + n_shift, ms_masses, tole_ptr);
      row_c[j] += matchMass(piece.suffix_mass[j] + c_shift, ms_masses, tole_ptr);
    }
  }
}

void addRows(std::vector<int>& row, const std::vector<int>& row_n,
             const std::vector<int>& row_c) {
  for (size_t j = 0; j < row.size(); j++) {
    row[j] += row_n[j] + row_c[j];
  }
}

// The site table of one PTM: s_table[i], for the residue i, is the number
// of fragments matched without the PTM before it plus those matched with
// the PTM after it (local_util::compOnePtmSTable, with the piece's known
// shifts and offsets).
std::vector<int> compOnePtmSTable(const Piece& piece,
                                  const ExtendMsPtrVec& extend_ms_ptr_vec,
                                  double ptm_mass, const LocalMngPtr& mng_ptr) {
  int len = piece.getLen();
  std::vector<int> row_1(len + 1, 0);
  std::vector<int> row_2(len + 1, 0);
  std::vector<int> row_n, row_c;
  // a site before the PTM: prefix without it, suffix with it
  compMatchRows(piece, extend_ms_ptr_vec, 0, ptm_mass, mng_ptr, row_n, row_c);
  addRows(row_1, row_n, row_c);
  // a site after the PTM: prefix with it, suffix without it
  compMatchRows(piece, extend_ms_ptr_vec, ptm_mass, 0, mng_ptr, row_n, row_c);
  addRows(row_2, row_n, row_c);
  std::vector<int> n_pref_score = local_util::compPrefScore(row_1);
  std::vector<int> c_suff_score = local_util::compSuffScore(row_2);
  std::vector<int> s_table(len, 0);
  for (int i = 0; i < len; i++) {
    s_table[i] = n_pref_score[i] + c_suff_score[i + 1];
  }
  return s_table;
}

// The three-layer match table of two PTMs (local_util::compTwoPtmMTable):
// m_table[0] before the first PTM, [1] between the PTMs, [2] after the
// second
std::vector<std::vector<int>> compTwoPtmMTable(
    const Piece& piece, const ExtendMsPtrVec& extend_ms_ptr_vec, double mass_1,
    double mass_2, const LocalMngPtr& mng_ptr) {
  int len = piece.getLen();
  std::vector<std::vector<int>> m_table(3, std::vector<int>(len + 1, 0));
  std::vector<int> row_n, row_c;
  compMatchRows(piece, extend_ms_ptr_vec, 0, mass_1 + mass_2, mng_ptr, row_n,
                row_c);
  addRows(m_table[0], row_n, row_c);
  compMatchRows(piece, extend_ms_ptr_vec, mass_1, mass_2, mng_ptr, row_n,
                row_c);
  addRows(m_table[1], row_n, row_c);
  compMatchRows(piece, extend_ms_ptr_vec, mass_1 + mass_2, 0, mng_ptr, row_n,
                row_c);
  addRows(m_table[2], row_n, row_c);
  return m_table;
}

// the best score of two PTMs in this order (LocalProcessor::compTwoPtmScr)
int compTwoPtmScr(const Piece& piece, const std::vector<std::vector<int>>& m_table,
                  const PtmPtr& ptm_ptr_1, const PtmPtr& ptm_ptr_2,
                  const Para& para) {
  int len = piece.getLen();
  std::vector<std::vector<int>> d_table(3, std::vector<int>(len + 1, 0));
  d_table[1][0] = NEG_SCORE;
  d_table[2][0] = NEG_SCORE;
  for (int i = 1; i <= len; i++) {
    d_table[0][i] = d_table[0][i - 1] + m_table[0][i];
    d_table[1][i] = d_table[1][i - 1] + m_table[1][i];
    if (modifiable(piece, i - 1, ptm_ptr_1, para) &&
        d_table[0][i - 1] > d_table[1][i - 1]) {
      d_table[1][i] = d_table[0][i - 1] + m_table[1][i];
    }
    d_table[2][i] = d_table[2][i - 1] + m_table[2][i];
    if (modifiable(piece, i - 1, ptm_ptr_2, para) &&
        d_table[1][i - 1] > d_table[2][i - 1]) {
      d_table[2][i] = d_table[1][i - 1] + m_table[2][i];
    }
  }
  return d_table[2][len];
}

bool hasSite(const Piece& piece, const PtmPtr& ptm_ptr, const Para& para) {
  for (int i = 0; i < piece.getLen(); i++) {
    if (modifiable(piece, i, ptm_ptr, para)) return true;
  }
  return false;
}

// a site of the first PTM before a site of the second
bool hasTwoSites(const Piece& piece, const PtmPtr& ptm_ptr_1,
                 const PtmPtr& ptm_ptr_2, const Para& para) {
  int len = piece.getLen();
  int first = -1;
  for (int i = 0; i < len; i++) {
    if (modifiable(piece, i, ptm_ptr_1, para)) {
      first = i;
      break;
    }
  }
  if (first < 0) return false;
  for (int i = first + 1; i < len; i++) {
    if (modifiable(piece, i, ptm_ptr_2, para)) return true;
  }
  return false;
}

bool massMatch(double mass, double shift_mass, double err_tole,
               const Para& para) {
  double diff = std::abs(mass - shift_mass);
  if (diff < err_tole) return true;
  return para.iso &&
         std::abs(diff - mass_constant::getIsotopeMass()) < err_tole;
}

// The site scores cut to the reported range (local_util::scrFilter with
// the threshold: the contiguous sites above it around the best site, the
// sum of their scores as the confidence); when no site reaches the
// threshold, the range of all the sites with a positive score, with the
// best site's score as the confidence (their sum would be 1.0 whatever
// the spread)
std::vector<double> filterScores(const std::vector<double>& scr_vec,
                                 double threshold, int& bgn, int& end,
                                 double& conf) {
  std::vector<double> scr = scr_vec;
  local_util::scrFilter(scr, bgn, end, conf, threshold);
  if (bgn >= 0) return scr;
  bgn = -1;
  end = -1;
  conf = 0;
  for (size_t i = 0; i < scr_vec.size(); i++) {
    if (scr_vec[i] > 0) {
      if (bgn < 0) bgn = i;
      end = i;
      conf = std::max(conf, scr_vec[i]);
    }
  }
  return std::vector<double>(scr_vec.begin() + bgn, scr_vec.begin() + end + 1);
}

MassShiftPtr makeShift(const Piece& piece, const PtmPtr& ptm_ptr,
                       const std::vector<double>& scr_vec, int score,
                       double threshold) {
  int bgn, end;
  double conf;
  std::vector<double> scr = filterScores(scr_vec, threshold, bgn, end, conf);
  LocalAnnoPtr anno = std::make_shared<LocalAnno>(
      bgn + piece.bgn, end + piece.bgn, conf, scr, score, ptm_ptr);
  AminoAcidPtr empty_aa_ptr = AminoAcidBase::getEmptyAminoAcidPtr();
  ModPtr mod_ptr = std::make_shared<Mod>(
      ResidueBase::getEmptyResiduePtr(),
      ResidueBase::getBaseResiduePtr(empty_aa_ptr, ptm_ptr),
      ModType::SIDE_CHAIN);
  AlterPtr alter = std::make_shared<Alter>(
      anno->getLeftBpPos(), anno->getRightBpPos() + 1, AlterType::VARIABLE,
      ptm_ptr->getMonoMass(), mod_ptr);
  alter->setLocalAnno(anno);
  return std::make_shared<MassShift>(alter);
}

// pow(prob_ratio, score - max_score): the normalization removes the
// factor, and the largest exponent stays 0
std::vector<double> compSiteProbs(const std::vector<int>& s_table,
                                  const std::vector<bool>& site_ok,
                                  double prob_ratio, int max_score) {
  std::vector<double> scr_vec(s_table.size(), 0.0);
  for (size_t i = 0; i < s_table.size(); i++) {
    if (site_ok[i]) {
      scr_vec[i] = std::pow(prob_ratio, s_table[i] - max_score);
    }
  }
  local_util::normalize(scr_vec);
  return scr_vec;
}

MassShiftPtrVec localizeOnePtm(const Piece& piece, const PtmPtr& ptm_ptr,
                               const ExtendMsPtrVec& extend_ms_ptr_vec,
                               const Para& para, const LocalMngPtr& mng_ptr,
                               int& best_score) {
  int len = piece.getLen();
  std::vector<int> s_table =
      compOnePtmSTable(piece, extend_ms_ptr_vec, ptm_ptr->getMonoMass(), mng_ptr);
  std::vector<bool> site_ok(len, false);
  best_score = 0;
  for (int i = 0; i < len; i++) {
    site_ok[i] = modifiable(piece, i, ptm_ptr, para);
    if (site_ok[i] && s_table[i] > best_score) best_score = s_table[i];
  }
  std::vector<double> scr_vec =
      compSiteProbs(s_table, site_ok, mng_ptr->prob_ratio_, best_score);
  MassShiftPtrVec result;
  result.push_back(
      makeShift(piece, ptm_ptr, scr_vec, best_score, mng_ptr->threshold_));
  return result;
}

// LocalProcessor::twoPtmLocalize on the piece: the probabilities of the
// site pairs, the break point of the highest probability of a pair around
// it, and the two site vectors on its sides
MassShiftPtrVec localizeTwoPtms(const Piece& piece, const PtmPtr& ptm_ptr_1,
                                const PtmPtr& ptm_ptr_2,
                                const std::vector<std::vector<int>>& m_table,
                                const Para& para, const LocalMngPtr& mng_ptr,
                                int& best_score) {
  int len = piece.getLen();
  std::vector<int> row_0 = m_table[0];
  std::vector<int> row_2 = m_table[2];
  std::vector<int> n_term_scores = local_util::compPrefScore(row_0);
  std::vector<int> c_term_scores = local_util::compSuffScore(row_2);
  // mid_pref[j]: the matches of the middle layer up to the site j
  std::vector<int> mid_pref(len + 1, 0);
  for (int j = 0; j <= len; j++) {
    mid_pref[j] = (j == 0 ? 0 : mid_pref[j - 1]) + m_table[1][j];
  }
  std::vector<int> pos_list_1;
  std::vector<int> pos_list_2;
  for (int i = 0; i < len; i++) {
    if (modifiable(piece, i, ptm_ptr_1, para)) pos_list_1.push_back(i);
    if (modifiable(piece, i, ptm_ptr_2, para)) pos_list_2.push_back(i);
  }
  // the score of the first PTM at pos_1 and the second at pos_2 > pos_1:
  // the N-terminal layer up to pos_1, the middle layer over (pos_1, pos_2],
  // the C-terminal layer from pos_2 + 1
  auto pairScore = [&](int pos_1, int pos_2) {
    return n_term_scores[pos_1] + c_term_scores[pos_2 + 1] +
           (mid_pref[pos_2] - mid_pref[pos_1]);
  };
  best_score = NEG_SCORE;
  for (int pos_1 : pos_list_1) {
    for (int pos_2 : pos_list_2) {
      if (pos_2 > pos_1) best_score = std::max(best_score, pairScore(pos_1, pos_2));
    }
  }
  std::vector<std::vector<double>> b_table(len, std::vector<double>(len, 0.0));
  double total_prob = 0.0;
  for (int pos_1 : pos_list_1) {
    for (int pos_2 : pos_list_2) {
      if (pos_2 > pos_1) {
        b_table[pos_1][pos_2] =
            std::pow(mng_ptr->prob_ratio_, pairScore(pos_1, pos_2) - best_score);
        total_prob += b_table[pos_1][pos_2];
      }
    }
  }
  if (total_prob <= 0.0) {
    LOG_ERROR("No valid modification sites!");
    return MassShiftPtrVec();
  }
  for (int pos_1 : pos_list_1) {
    for (int pos_2 : pos_list_2) {
      if (pos_2 > pos_1) b_table[pos_1][pos_2] /= total_prob;
    }
  }
  // the break point: the probability that the first PTM is at or before it
  // and the second after it
  std::vector<double> break_scores(len - 1, 0.0);
  for (int pos_1 : pos_list_1) {
    for (int pos_2 : pos_list_2) {
      if (pos_2 > pos_1) {
        for (int p = pos_1; p < pos_2 && p < len - 1; p++) {
          break_scores[p] += b_table[pos_1][pos_2];
        }
      }
    }
  }
  int break_pos = std::distance(
      break_scores.begin(),
      std::max_element(break_scores.begin(), break_scores.end()));
  std::vector<double> scr_vec_1(len, 0.0);
  std::vector<double> scr_vec_2(len, 0.0);
  for (int pos_1 : pos_list_1) {
    if (pos_1 > break_pos) continue;
    for (int pos_2 : pos_list_2) {
      if (pos_2 > break_pos) {
        scr_vec_1[pos_1] += b_table[pos_1][pos_2];
        scr_vec_2[pos_2] += b_table[pos_1][pos_2];
      }
    }
  }
  MassShiftPtrVec result;
  result.push_back(
      makeShift(piece, ptm_ptr_1, scr_vec_1, best_score, mng_ptr->threshold_));
  result.push_back(
      makeShift(piece, ptm_ptr_2, scr_vec_2, best_score, mng_ptr->threshold_));
  return result;
}

}  // namespace

std::string Explanation::getName() const {
  std::string name = first->getAbbrName();
  if (second != nullptr) {
    name += "+" + second->getAbbrName();
  }
  return name;
}

bool Explanation::hasPtm(const std::string& abbr) const {
  if (first->getAbbrName() == abbr) return true;
  return second != nullptr && second->getAbbrName() == abbr;
}

std::vector<size_t> breakTies(const std::vector<Explanation>& expls,
                              const std::vector<bool>& acetyl_forms,
                              std::string& rule) {
  rule = "";
  std::vector<size_t> kept;
  for (size_t i = 0; i < expls.size(); i++) kept.push_back(i);
  auto apply = [&](const std::string& name,
                   const std::function<bool(size_t)>& drop) {
    if (kept.size() <= 1) return;
    std::vector<size_t> next;
    for (size_t k = 0; k < kept.size(); k++) {
      if (!drop(kept[k])) next.push_back(kept[k]);
    }
    if (next.size() < kept.size() && !next.empty()) {
      kept = next;
      rule += (rule.empty() ? "" : "+") + name;
    }
  };
  // (1) fewer PTMs
  int min_num = 2;
  for (size_t k = 0; k < kept.size(); k++) {
    min_num = std::min(min_num, expls[kept[k]].getPtmNum());
  }
  apply("fewer_ptms",
        [&](size_t i) { return expls[i].getPtmNum() > min_num; });
  // (2) the N-terminal acetylation form over an Acetyl on the plain form
  bool any_acetyl_form = false;
  for (size_t k = 0; k < kept.size(); k++) {
    if (acetyl_forms[kept[k]]) any_acetyl_form = true;
  }
  if (any_acetyl_form) {
    apply("acetyl_form", [&](size_t i) {
      return !acetyl_forms[i] && expls[i].hasPtm("Acetyl");
    });
  }
  return kept;
}

bool isKnownShift(const MassShiftPtr& shift_ptr) {
  AlterTypePtr type_ptr = shift_ptr->getTypePtr();
  return type_ptr == AlterType::VARIABLE || type_ptr == AlterType::UNEXPECTED;
}

Piece makePiece(const ProteoformPtr& form_ptr, int bgn, int end,
                const MassShiftPtr& shift_ptr) {
  Piece piece;
  piece.bgn = bgn;
  piece.end = end;
  piece.form_len = form_ptr->getLen();
  const ResiduePtrVec& residues = form_ptr->getResSeqPtr()->getResidues();
  piece.residues.assign(residues.begin() + bgn, residues.begin() + end);
  if (shift_ptr != nullptr) {
    piece.shift_mass = shift_ptr->getMassShift();
  }
  for (int i = 0; i < bgn; i++) {
    piece.n_offset += residues[i]->getMass();
  }
  for (int i = end; i < piece.form_len; i++) {
    piece.c_offset += residues[i]->getMass();
  }
  const MassShiftPtrVec& shifts = form_ptr->getMassShiftPtrVec();
  for (size_t k = 0; k < shifts.size(); k++) {
    MassShiftPtr s = shifts[k];
    if (s == shift_ptr) continue;
    int left = s->getLeftBpPos();
    int right = s->getRightBpPos();
    if (s->getTypePtr() == AlterType::FIXED) {
      if (left >= bgn && right <= end) {
        piece.fixed_shifts.push_back({left - bgn, right - bgn, s->getMassShift()});
      }
      continue;
    }
    if (!isKnownShift(s)) continue;
    // a shift belongs to the piece of its left position; one that reaches
    // over the piece's end is cut at it
    if (left < bgn) {
      piece.n_offset += s->getMassShift();
    } else if (left >= end) {
      piece.c_offset += s->getMassShift();
    } else {
      piece.known_shifts.push_back(
          {left - bgn, std::min(right, end) - bgn, s->getMassShift()});
    }
  }
  int len = piece.getLen();
  piece.prefix_mass.assign(len + 1, 0.0);
  piece.suffix_mass.assign(len + 1, 0.0);
  piece.valid.assign(len + 1, true);
  // The protein N-terminal modification (the acetylation of a proteoform
  // read from the PrSM file) is the residue sequence's N-terminal mod, not
  // part of any residue's mass: it belongs to the prefix masses of a piece
  // that starts at the N-terminus and to the N offset of any other piece
  // (as BpSpec starts its prefix masses with it). Until this was added
  // (toppic_test experiment 140, amendment 7) the N-terminal fragments of
  // acetylated proteoforms were 42.01 Da off and never matched.
  double n_mod_mass = form_ptr->getResSeqPtr()->getNModPtr()->getShift();
  if (bgn > 0) {
    piece.n_offset += n_mod_mass;
    n_mod_mass = 0;
  }
  for (int j = 1; j <= len; j++) {
    piece.prefix_mass[j] =
        piece.prefix_mass[j - 1] + piece.residues[j - 1]->getMass();
  }
  for (int j = len - 1; j >= 0; j--) {
    piece.suffix_mass[j] = piece.suffix_mass[j + 1] + piece.residues[j]->getMass();
  }
  for (int j = 1; j <= len; j++) piece.prefix_mass[j] += n_mod_mass;
  piece.suffix_mass[0] += n_mod_mass;
  for (size_t k = 0; k < piece.known_shifts.size(); k++) {
    const ShiftRange& r = piece.known_shifts[k];
    for (int j = 0; j <= len; j++) {
      if (r.right <= j) piece.prefix_mass[j] += r.mass;
      if (r.left >= j) piece.suffix_mass[j] += r.mass;
      if (r.left < j && j < r.right) piece.valid[j] = false;
    }
  }
  return piece;
}

std::vector<int> compMatchedSites(const ProteoformPtr& form_ptr,
                                  const ExtendMsPtrVec& extend_ms_ptr_vec,
                                  const LocalMngPtr& mng_ptr) {
  // every shift of the proteoform is known here
  Piece piece = makePiece(form_ptr, 0, form_ptr->getLen(), nullptr);
  std::vector<int> row_n, row_c;
  compMatchRows(piece, extend_ms_ptr_vec, 0, 0, mng_ptr, row_n, row_c);
  std::vector<int> matched(piece.getLen() + 1, 0);
  for (size_t j = 0; j < matched.size(); j++) {
    matched[j] = (row_n[j] + row_c[j] > 0) ? 1 : 0;
  }
  return matched;
}

int countMatchedSites(const std::vector<int>& matched_sites, int from, int to) {
  int n = 0;
  for (int j = std::max(from, 0); j <= to && j < static_cast<int>(matched_sites.size()); j++) {
    n += matched_sites[j];
  }
  return n;
}

MassShiftPtr shiftPositions(const MassShiftPtr& shift_ptr, int delta) {
  if (delta == 0) return shift_ptr;
  AlterPtr ori = shift_ptr->getAlterPtr(0);
  AlterPtr alter = std::make_shared<Alter>(
      ori->getLeftBpPos() - delta, ori->getRightBpPos() - delta,
      ori->getTypePtr(), ori->getMass(), ori->getModPtr());
  LocalAnnoPtr anno = ori->getLocalAnno();
  if (anno != nullptr) {
    alter->setLocalAnno(std::make_shared<LocalAnno>(
        anno->getLeftBpPos() - delta, anno->getRightBpPos() - delta,
        anno->getConf(), anno->getScrVec(), anno->getRawScr(),
        anno->getPtmPtr()));
  }
  return std::make_shared<MassShift>(alter);
}

int findBreakSite(const MassShiftPtr& shift_1, const MassShiftPtr& shift_2,
                  const std::vector<int>& matched_sites, int& matched_num) {
  int lo = shift_1->getRightBpPos();
  int hi = shift_2->getLeftBpPos();
  std::vector<int> sites;
  for (int j = lo; j <= hi; j++) {
    if (matched_sites[j]) sites.push_back(j);
  }
  matched_num = sites.size();
  if (sites.empty()) {
    return (lo + hi) / 2;
  }
  return sites[sites.size() / 2];
}

bool modifiable(const Piece& piece, int i, const PtmPtr& ptm_ptr,
                const Para& para) {
  for (size_t k = 0; k < piece.fixed_shifts.size(); k++) {
    if (piece.fixed_shifts[k].left <= i && i < piece.fixed_shifts[k].right) {
      return false;
    }
  }
  const ResiduePtr& residue_ptr = piece.residues[i];
  // the N-term and C-term entries of the file apply to the proteoform's
  // first and last residue; the any entries are in all three lists
  int form_pos = piece.bgn + i;
  const ModPtrVec* mod_list = &para.mod_list_any;
  if (form_pos == 0) {
    mod_list = &para.mod_list_N;
  } else if (form_pos == piece.form_len - 1) {
    mod_list = &para.mod_list_C;
  }
  for (size_t j = 0; j < mod_list->size(); j++) {
    const ModPtr& mod_ptr = (*mod_list)[j];
    if (mod_ptr->getOriResiduePtr()->isSame(residue_ptr) &&
        mod_ptr->getModResiduePtr()->getPtmPtr()->isSame(ptm_ptr)) {
      return true;
    }
  }
  return false;
}

std::vector<Explanation> findExplanations(const Piece& piece, double err_tole,
                                          const Para& para) {
  std::vector<Explanation> result;
  const PtmPtrVec& ptms = para.ptm_ptr_vec;
  for (size_t i = 0; i < ptms.size(); i++) {
    if (massMatch(ptms[i]->getMonoMass(), piece.shift_mass, err_tole, para) &&
        hasSite(piece, ptms[i], para)) {
      result.push_back({ptms[i], nullptr});
    }
  }
  if (para.max_ptm < 2) return result;
  for (size_t i = 0; i < ptms.size(); i++) {
    for (size_t j = i; j < ptms.size(); j++) {
      double mass = ptms[i]->getMonoMass() + ptms[j]->getMonoMass();
      if (!massMatch(mass, piece.shift_mass, err_tole, para)) continue;
      if (hasTwoSites(piece, ptms[i], ptms[j], para) ||
          hasTwoSites(piece, ptms[j], ptms[i], para)) {
        result.push_back({ptms[i], ptms[j]});
      }
    }
  }
  return result;
}

MassShiftPtrVec localize(const Piece& piece, const Explanation& expl,
                         const ExtendMsPtrVec& extend_ms_ptr_vec,
                         const Para& para, const LocalMngPtr& mng_ptr,
                         int& best_score) {
  if (expl.second == nullptr) {
    return localizeOnePtm(piece, expl.first, extend_ms_ptr_vec, para, mng_ptr,
                          best_score);
  }
  // the order of the two PTMs with the higher score
  PtmPtr ptm_1 = expl.first;
  PtmPtr ptm_2 = expl.second;
  std::vector<std::vector<int>> m_table = compTwoPtmMTable(
      piece, extend_ms_ptr_vec, ptm_1->getMonoMass(), ptm_2->getMonoMass(),
      mng_ptr);
  int score = NEG_SCORE;
  if (hasTwoSites(piece, ptm_1, ptm_2, para)) {
    score = compTwoPtmScr(piece, m_table, ptm_1, ptm_2, para);
  }
  if (!ptm_1->isSame(ptm_2) && hasTwoSites(piece, ptm_2, ptm_1, para)) {
    std::vector<std::vector<int>> m_table_2 = compTwoPtmMTable(
        piece, extend_ms_ptr_vec, ptm_2->getMonoMass(), ptm_1->getMonoMass(),
        mng_ptr);
    int score_2 = compTwoPtmScr(piece, m_table_2, ptm_2, ptm_1, para);
    if (score_2 > score) {
      std::swap(ptm_1, ptm_2);
      m_table = m_table_2;
      score = score_2;
    }
  }
  return localizeTwoPtms(piece, ptm_1, ptm_2, m_table, para, mng_ptr,
                         best_score);
}

}  // namespace local_align

}  // namespace toppic
