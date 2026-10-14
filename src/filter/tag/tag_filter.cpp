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

#include "filter/tag/tag_filter.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <set>
#include <tuple>
#include <unordered_map>
#include <utility>

#include "common/base/mass_constant.hpp"
#include "common/base/mod_util.hpp"
#include "common/base/prot_mod_util.hpp"
#include "common/base/proteoform_type.hpp"
#include "common/base/residue_base.hpp"
#include "common/util/logger.hpp"
#include "ms/factory/prm_ms_util.hpp"
#include "ms/spec/base_peak_type.hpp"
#include "seq/proteoform_factory.hpp"

namespace toppic {

namespace {

// the standard amino acids a tag can hold; I is read as L
const std::string kTagLetters = "ACDEFGHKLMNPQRSTVWY";

// The sums of up to num masses of the list (with repetition), 0 included,
// each sum once (VarPtmFilterMng::computeShifts, in 1e-4 Da units).
std::vector<double> computeShiftSums(const std::vector<double>& masses,
                                     int num) {
  const double scale = 10000;
  std::vector<long> round_sums = {0};
  for (int i = 0; i < num; i++) {
    size_t cur_len = round_sums.size();
    for (size_t j = 0; j < cur_len; j++) {
      for (size_t k = 0; k < masses.size(); k++) {
        long new_sum = round_sums[j] + std::lround(masses[k] * scale);
        if (std::find(round_sums.begin(), round_sums.end(), new_sum) ==
            round_sums.end()) {
          round_sums.push_back(new_sum);
        }
      }
    }
  }
  std::vector<double> sums;
  for (size_t i = 0; i < round_sums.size(); i++) {
    sums.push_back(round_sums[i] / scale);
  }
  return sums;
}

// the indexes i in [first, last] with |masses[i] - target| <= tole;
// masses is sorted
void findMasses(const std::vector<double>& masses, int first, int last,
                double target, double tole, std::vector<int>* results) {
  if (first > last) {
    return;
  }
  auto it = std::lower_bound(masses.begin() + first, masses.begin() + last + 1,
                             target - tole);
  while (it != masses.begin() + last + 1 && *it <= target + tole) {
    results->push_back(static_cast<int>(it - masses.begin()));
    ++it;
  }
}

}  // namespace

TagFilter::TagFilter(const TagIndexPtr& index_ptr,
                     const TagFilterMngPtr& mng_ptr)
    : index_ptr_(index_ptr), mng_ptr_(mng_ptr) {
  initResidues();
  if (mng_ptr_->stage_ == TagFilterMng::Stage::kVar) {
    initVarPtmShifts();
  }
}

// The residue masses the tags are read with: the 19 letters (I as L) with
// the fixed modifications applied (C57: C is 160.03 Da).
void TagFilter::initResidues() {
  ResiduePtrVec base_residues = ResidueBase::getBaseNonePtmResiduePtrVec();
  ResiduePtrVec residues = mod_util::geneResidueListWithMod(
      base_residues, mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
  for (size_t i = 0; i < residues.size(); i++) {
    const std::string& letter =
        residues[i]->getAminoAcidPtr()->getOneLetter();
    if (letter.size() != 1 || kTagLetters.find(letter[0]) == std::string::npos) {
      continue;
    }
    TagResidue res;
    res.letter = letter[0];
    res.mass = residues[i]->getMass();
    residues_.push_back(res);
  }
  std::sort(residues_.begin(), residues_.end(),
            [](const TagResidue& a, const TagResidue& b) {
              return a.mass < b.mass;
            });
  if (residues_.empty()) {
    LOG_ERROR("The tag filter found no residue masses!");
    exit(EXIT_FAILURE);
  }
  min_res_mass_ = residues_.front().mass;
  max_res_mass_ = residues_.back().mass;
}

void TagFilter::initVarPtmShifts() {
  std::vector<double> single_shifts;
  if (mng_ptr_->var_ptm_file_name_ != "") {
    single_shifts =
        mod_util::readModTxtToShiftList(mng_ptr_->var_ptm_file_name_);
    ModPtrVec n_term_mod_ptrs =
        mod_util::readNTermModTxt(mng_ptr_->var_ptm_file_name_);
    for (size_t i = 0; i < n_term_mod_ptrs.size(); i++) {
      double shift = n_term_mod_ptrs[i]->getShift();
      if (std::find(n_term_mod_shifts_.begin(), n_term_mod_shifts_.end(),
                    shift) == n_term_mod_shifts_.end()) {
        n_term_mod_shifts_.push_back(shift);
      }
    }
  }
  int num = std::max(mng_ptr_->var_ptm_num_, 0);
  var_shift_sums_ = computeShiftSums(single_shifts, num);
  var_shift_sums_fewer_ = computeShiftSums(single_shifts, std::max(num - 1, 0));
}

std::vector<SimplePrsmPtrVec> TagFilter::process(
    const SpectrumSetPtr& spec_set_ptr) const {
  Stats stats;
  return process(spec_set_ptr, &stats);
}

std::vector<SimplePrsmPtrVec> TagFilter::process(
    const SpectrumSetPtr& spec_set_ptr, Stats* stats) const {
  std::vector<SimplePrsmPtrVec> results(4);
  const PrmMsPtrVec& prm_ms_ptr_vec = spec_set_ptr->getMsTwoPtrVec();
  if (prm_ms_ptr_vec.empty()) {
    return results;
  }
  PeakTolerancePtr tole_ptr =
      mng_ptr_->prsm_para_ptr_->getSpParaPtr()->getPrecTolerancePtr();
  // the PRM peaks of the group's spectra with the zero and precursor
  // peaks, sorted, with their ids
  PrmPeakPtrVec peaks = prm_ms_util::getPrmPeakPtrs(prm_ms_ptr_vec, tole_ptr);
  std::vector<Candidate> cands = findCandidates(peaks, stats);
  // the header of the PRM spectra holds the label-free precursor mass, as
  // the mass filters' candidates do
  MsHeaderPtr header_ptr = prm_ms_ptr_vec[0]->getMsHeaderPtr();
  int group_spec_num = static_cast<int>(prm_ms_ptr_vec.size());
  for (size_t i = 0; i < cands.size(); i++) {
    switch (mng_ptr_->stage_) {
      case TagFilterMng::Stage::kOne:
        geneOneShiftPrsms(cands[i], header_ptr, group_spec_num, &results);
        break;
      case TagFilterMng::Stage::kVar:
        geneVarPtmPrsms(cands[i], header_ptr, group_spec_num, &results);
        break;
      case TagFilterMng::Stage::kMulti:
      case TagFilterMng::Stage::kTopmg:
      default:
        geneProtPrsm(cands[i], header_ptr, group_spec_num, &results);
        break;
    }
  }
  return results;
}

// Depth-first enumeration of the paths of tag_len_ edges from the current
// path (str, peak_ids; the last peak id is the path's end).
void TagFilter::extendTag(
    const std::vector<std::vector<std::pair<int, char>>>& edges,
    std::string* str, std::vector<int>* peak_ids,
    std::vector<Tag>* tags) const {
  if (static_cast<int>(str->size()) == mng_ptr_->tag_len_) {
    Tag tag;
    tag.str = *str;
    tag.peak_ids = *peak_ids;
    tags->push_back(std::move(tag));
    return;
  }
  // the group index of the path's last peak
  int last = peak_ids->back();
  const std::vector<std::pair<int, char>>& out = edges[last];
  for (size_t i = 0; i < out.size(); i++) {
    str->push_back(out[i].second);
    peak_ids->push_back(out[i].first);
    extendTag(edges, str, peak_ids, tags);
    peak_ids->pop_back();
    str->pop_back();
  }
}

void TagFilter::geneTags(const PrmPeakPtrVec& peaks,
                         std::vector<Tag>* tags) const {
  // the peaks of each base type, in increasing position (peaks is sorted)
  std::vector<int> original;
  std::vector<int> reversed;
  for (size_t i = 0; i < peaks.size(); i++) {
    if (peaks[i]->getBaseTypePtr() == BasePeakType::REVERSED) {
      reversed.push_back(static_cast<int>(i));
    } else {
      original.push_back(static_cast<int>(i));
    }
  }
  // peaks of one type within their tolerances of each other are one peak
  // (a fragment mass at several charge states); the first is kept
  auto dedup = [&](const std::vector<int>& all) {
    std::vector<int> kept;
    for (size_t i = 0; i < all.size(); i++) {
      const PrmPeakPtr& p = peaks[all[i]];
      if (!kept.empty()) {
        const PrmPeakPtr& prev = peaks[kept.back()];
        if (p->getPosition() - prev->getPosition() <=
            p->getStrictTolerance() + prev->getStrictTolerance()) {
          continue;
        }
      }
      kept.push_back(all[i]);
    }
    return kept;
  };
  std::vector<int> original_kept = dedup(original);
  std::vector<int> reversed_kept = dedup(reversed);
  const std::vector<int>* groups[2] = {&original_kept, &reversed_kept};
  for (int g = 0; g < 2; g++) {
    const std::vector<int>& group = *groups[g];
    size_t n = group.size();
    // edges[a]: the peaks b > a of the group one residue above peak a, with
    // the residue letter; the edge targets are peak indexes of peaks
    std::vector<std::vector<std::pair<int, char>>> edges(peaks.size());
    for (size_t a = 0; a < n; a++) {
      const PrmPeakPtr& pa = peaks[group[a]];
      double pos_a = pa->getPosition();
      double tole_a = pa->getStrictTolerance();
      for (size_t b = a + 1; b < n; b++) {
        const PrmPeakPtr& pb = peaks[group[b]];
        double diff = pb->getPosition() - pos_a;
        double tole = tole_a + pb->getStrictTolerance();
        if (diff > max_res_mass_ + tole) {
          break;
        }
        if (diff < min_res_mass_ - tole) {
          continue;
        }
        for (size_t r = 0; r < residues_.size(); r++) {
          if (std::abs(diff - residues_[r].mass) <= tole) {
            edges[group[a]].push_back(
                std::make_pair(group[b], residues_[r].letter));
          }
        }
      }
    }
    for (size_t a = 0; a < n; a++) {
      if (edges[group[a]].empty()) {
        continue;
      }
      std::string str;
      std::vector<int> peak_ids = {group[a]};
      extendTag(edges, &str, &peak_ids, tags);
    }
  }
  // the peak indexes are the peak ids (getPrmPeakPtrs numbers the sorted
  // peaks), so nothing to convert
}

// The diagonals of a protein: its occurrences sorted by shift, a new
// diagonal where two neighbours are further apart than the sum of their
// first peaks' tolerances; a diagonal's support is the number of distinct
// peaks its tags cover. The score is the support of the best diagonal plus
// that of the second when at least two tags support it
// (TagFilterMng::getMinScore), so that two chance tags on different
// diagonals of a long protein do not add up.
void TagFilter::scoreHit(ProtHit* hit) const {
  std::vector<Occurrence>& occs = hit->occs;
  std::sort(occs.begin(), occs.end(),
            [](const Occurrence& a, const Occurrence& b) {
              return a.shift < b.shift;
            });
  hit->diagonals.clear();
  size_t i = 0;
  while (i < occs.size()) {
    size_t j = i + 1;
    while (j < occs.size() &&
           occs[j].shift - occs[j - 1].shift <=
               occs[j].first_tole + occs[j - 1].first_tole) {
      j++;
    }
    double sum = 0;
    std::set<int> peak_ids;
    for (size_t k = i; k < j; k++) {
      sum += occs[k].shift;
      peak_ids.insert(occs[k].peak_ids.begin(), occs[k].peak_ids.end());
    }
    Diagonal diagonal;
    diagonal.shift = sum / static_cast<double>(j - i);
    diagonal.support = static_cast<int>(peak_ids.size());
    hit->diagonals.push_back(diagonal);
    i = j;
  }
  std::sort(hit->diagonals.begin(), hit->diagonals.end(),
            [](const Diagonal& a, const Diagonal& b) {
              if (a.support != b.support) {
                return a.support > b.support;
              }
              return std::abs(a.shift) < std::abs(b.shift);
            });
  // the best diagonal, and the second when two or more tags support it
  hit->score = 0;
  if (!hit->diagonals.empty()) {
    hit->score = hit->diagonals[0].support;
    if (hit->diagonals.size() > 1 &&
        hit->diagonals[1].support >= mng_ptr_->getMinScore()) {
      hit->score += hit->diagonals[1].support;
    }
  }
}

std::vector<TagFilter::Candidate> TagFilter::findCandidates(
    const PrmPeakPtrVec& peaks, Stats* stats) const {
  std::vector<Candidate> cands;
  std::vector<Tag> tags;
  geneTags(peaks, &tags);
  stats->tag_num = tags.size();
  if (tags.empty()) {
    return cands;
  }
  // one lookup per distinct tag string
  std::sort(tags.begin(), tags.end(),
            [](const Tag& a, const Tag& b) { return a.str < b.str; });
  std::unordered_map<int, ProtHit> hits;
  std::vector<TagIndex::Match> matches;
  size_t occurrence_num = 0;
  size_t i = 0;
  while (i < tags.size()) {
    size_t j = i;
    while (j < tags.size() && tags[j].str == tags[i].str) {
      j++;
    }
    index_ptr_->findTag(tags[i].str, &matches);
    occurrence_num += matches.size() * (j - i);
    for (size_t m = 0; m < matches.size(); m++) {
      ProtHit& hit = hits[matches[m].prot_idx];
      double prefix_mass =
          index_ptr_->getPrefixMass(matches[m].prot_idx, matches[m].pos);
      for (size_t t = i; t < j; t++) {
        const Tag& tag = tags[t];
        Occurrence occ;
        occ.pos = matches[m].pos;
        const PrmPeakPtr& first = peaks[tag.peak_ids[0]];
        occ.first_mass = first->getPosition();
        occ.first_tole = first->getNStrictCRelaxTolerance();
        occ.shift = occ.first_mass - prefix_mass;
        occ.peak_ids = tag.peak_ids;
        hit.occs.push_back(std::move(occ));
      }
    }
    i = j;
  }
  stats->occurrence_num = occurrence_num;
  stats->prot_num = hits.size();
  // the proteins by score, then name
  std::vector<std::pair<int, int>> order;  // (prot_idx, score)
  int min_score = mng_ptr_->getMinScore();
  for (auto it = hits.begin(); it != hits.end(); ++it) {
    scoreHit(&it->second);
    if (it->second.score >= min_score) {
      order.push_back(std::make_pair(it->first, it->second.score));
    }
  }
  // by score, then the shorter protein (the same tags are less likely by
  // chance in a shorter protein), then the name
  std::sort(order.begin(), order.end(),
            [this](const std::pair<int, int>& a,
                   const std::pair<int, int>& b) {
              if (a.second != b.second) {
                return a.second > b.second;
              }
              int len_a = index_ptr_->getProtLen(a.first);
              int len_b = index_ptr_->getProtLen(b.first);
              if (len_a != len_b) {
                return len_a < len_b;
              }
              return index_ptr_->getFastaSeqPtr(a.first)->getName() <
                     index_ptr_->getFastaSeqPtr(b.first)->getName();
            });
  for (size_t k = 0; k < order.size() &&
                     static_cast<int>(k) < mng_ptr_->cand_num_;
       k++) {
    Candidate cand;
    cand.prot_idx = order[k].first;
    cand.score = order[k].second;
    ProtHit& hit = hits[order[k].first];
    cand.occs = std::move(hit.occs);
    cand.diagonals = std::move(hit.diagonals);
    cands.push_back(std::move(cand));
  }
  return cands;
}

// one stage: the tag diagonals as the candidate's shifts, the same
// candidate in every type
void TagFilter::geneOneShiftPrsms(const Candidate& cand,
                                  const MsHeaderPtr& header_ptr,
                                  int group_spec_num,
                                  std::vector<SimplePrsmPtrVec>* results) const {
  const FastaSeqPtr& seq_ptr = index_ptr_->getFastaSeqPtr(cand.prot_idx);
  ProteoformPtr proteo_ptr = proteoform_factory::geneDbProteoformPtr(
      seq_ptr, mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
  // the best supported diagonals; two diagonals closer than
  // shift_merge_tole_ are one shift for the search
  double prec_mass = header_ptr->getFirstPrecMonoMass();
  double prot_mass = proteo_ptr->getResSeqPtr()->getSeqMass();
  std::vector<double> n_term_shifts;
  std::vector<double> c_term_shifts;
  for (size_t k = 0; k < cand.diagonals.size() &&
                     static_cast<int>(n_term_shifts.size()) <
                         mng_ptr_->max_shift_num_;
       k++) {
    double shift = cand.diagonals[k].shift;
    bool listed = false;
    for (size_t j = 0; j < n_term_shifts.size() && !listed; j++) {
      listed = std::abs(n_term_shifts[j] - shift) <= mng_ptr_->shift_merge_tole_;
    }
    if (listed) {
      continue;
    }
    n_term_shifts.push_back(shift);
    // SimplePrsm::getNTermShiftsFromCTermShifts gives prec - (prot + c)
    c_term_shifts.push_back(prec_mass - prot_mass - shift);
  }
  SimplePrsmPtr prsm_ptr = std::make_shared<SimplePrsm>(
      header_ptr, group_spec_num, proteo_ptr, cand.score);
  prsm_ptr->setNTermShifts(n_term_shifts);
  prsm_ptr->setCTermShifts(c_term_shifts);
  const PrsmParaPtr& prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  const std::vector<ProteoformTypePtr> types = {
      ProteoformType::COMPLETE, ProteoformType::PREFIX, ProteoformType::SUFFIX,
      ProteoformType::INTERNAL};
  for (size_t t = 0; t < types.size(); t++) {
    if (prsm_para_ptr->allowProtType(types[t])) {
      (*results)[t].push_back(prsm_ptr);
    }
  }
}

// var stage: the proteoforms of the protein that contain a tag and whose
// mass matches the precursor up to the variable PTMs
void TagFilter::geneVarPtmPrsms(const Candidate& cand,
                                const MsHeaderPtr& header_ptr,
                                int group_spec_num,
                                std::vector<SimplePrsmPtrVec>* results) const {
  const PrsmParaPtr& prsm_para_ptr = mng_ptr_->prsm_para_ptr_;
  const FastaSeqPtr& seq_ptr = index_ptr_->getFastaSeqPtr(cand.prot_idx);
  ProteoformPtr proteo_ptr = proteoform_factory::geneDbProteoformPtr(
      seq_ptr, prsm_para_ptr->getFixModPtrVec());
  std::vector<double> prms = proteo_ptr->getBpSpecPtr()->getPrmMasses();
  const ResiduePtrVec& residues = proteo_ptr->getResSeqPtr()->getResidues();
  int len = static_cast<int>(residues.size());
  if (static_cast<int>(prms.size()) != len + 1) {
    LOG_ERROR("The prefix masses of " << seq_ptr->getName()
                                      << " do not match its residues!");
    return;
  }
  double prec_mass = header_ptr->getFirstPrecMonoMass();
  double prec_tole =
      prsm_para_ptr->getSpParaPtr()->getPrecTolerancePtr()->compPrecErrorTole(
          prec_mass);
  double water = mass_constant::getWaterMass();
  const ProtModPtrVec& prot_mods = prsm_para_ptr->getProtModPtrVec();
  // the forms the protein allows, checked once
  std::vector<bool> form_allowed(prot_mods.size());
  for (size_t f = 0; f < prot_mods.size(); f++) {
    form_allowed[f] = prot_mod_util::allowMod(prot_mods[f], residues);
  }
  // (start, end, form index or -1, N-term mod index or -1)
  std::set<std::tuple<int, int, int, int>> seen;
  int tag_len = mng_ptr_->tag_len_;
  std::vector<int> starts;
  std::vector<int> ends;
  // the proteoform [start, end] with the N-terminal masses n_extra (the
  // form's modification mass or the truncation N-terminal modification)
  auto addProteoform = [&](int start, int end, int form_idx, int mod_idx,
                           double n_term_shift, double trunc_shift) {
    std::tuple<int, int, int, int> key(start, end, form_idx, mod_idx);
    if (!seen.insert(key).second) {
      return;
    }
    ProteoformTypePtr type_ptr;
    if (form_idx >= 0) {
      type_ptr = end == len - 1 ? ProteoformType::COMPLETE
                                : ProteoformType::PREFIX;
    } else {
      type_ptr = end == len - 1 ? ProteoformType::SUFFIX
                                : ProteoformType::INTERNAL;
    }
    if (!prsm_para_ptr->allowProtType(type_ptr)) {
      return;
    }
    SimplePrsmPtr prsm_ptr = std::make_shared<SimplePrsm>(
        header_ptr, group_spec_num, proteo_ptr, cand.score);
    prsm_ptr->setSeqRange(start, end);
    prsm_ptr->setNTermShifts(std::vector<double>(1, n_term_shift));
    prsm_ptr->setCTermShifts(std::vector<double>(1, 0.0));
    if (mod_idx >= 0) {
      prsm_ptr->setTruncNTermModShifts(std::vector<double>(1, trunc_shift));
    }
    int t = type_ptr == ProteoformType::COMPLETE ? 0
            : type_ptr == ProteoformType::PREFIX  ? 1
            : type_ptr == ProteoformType::SUFFIX  ? 2
                                                  : 3;
    (*results)[t].push_back(prsm_ptr);
  };
  for (size_t i = 0; i < cand.occs.size(); i++) {
    const Occurrence& occ = cand.occs[i];
    if (occ.pos < 0 || occ.pos + tag_len > len) {
      continue;
    }
    double s = occ.shift;
    double tole = occ.first_tole;
    // the proteoform's last residue is at or after the tag's last one
    int first_end_bp = occ.pos + tag_len;
    // (a) a protein N-terminal form: start at the form's position, the
    // N-terminal diagonal at the form's shift plus variable PTMs before
    // the tag
    for (size_t f = 0; f < prot_mods.size(); f++) {
      if (!form_allowed[f]) {
        continue;
      }
      int start = prot_mods[f]->getTruncPtr()->getTruncLen();
      if (start > occ.pos) {
        continue;
      }
      double form_shift = prot_mods[f]->getProtShift();
      bool fits = false;
      for (size_t c = 0; c < var_shift_sums_.size() && !fits; c++) {
        fits = std::abs(s - (form_shift + var_shift_sums_[c])) <= tole;
      }
      if (!fits) {
        continue;
      }
      // the mass of the form's modification (the acetylation; 0 for NONE
      // and NME, whose removed M is the start)
      double pep_shift = prot_mods[f]->getPepShift();
      for (size_t c = 0; c < var_shift_sums_.size(); c++) {
        double target = prec_mass - var_shift_sums_[c] - water - pep_shift +
                        prms[start];
        ends.clear();
        findMasses(prms, first_end_bp, len, target, prec_tole, &ends);
        for (size_t e = 0; e < ends.size(); e++) {
          addProteoform(start, ends[e] - 1, static_cast<int>(f), -1,
                        form_shift, 0.0);
        }
      }
    }
    // (b) a truncation start i: -prms[i] + variable PTMs before the tag
    // (+ a truncation N-terminal modification) = s
    for (int m = -1; m < static_cast<int>(n_term_mod_shifts_.size()); m++) {
      double trunc_shift = m < 0 ? 0.0 : n_term_mod_shifts_[m];
      const std::vector<double>& sums =
          m < 0 ? var_shift_sums_ : var_shift_sums_fewer_;
      // a truncation N-terminal modification needs a start at residue 3
      // or later (Proteoform::setTruncNTermModPtr)
      int first_start = m < 0 ? 1 : 2;
      for (size_t c = 0; c < sums.size(); c++) {
        double target_start = sums[c] + trunc_shift - s;
        starts.clear();
        findMasses(prms, first_start, occ.pos, target_start, tole, &starts);
        for (size_t k = 0; k < starts.size(); k++) {
          int start = starts[k];
          for (size_t c2 = 0; c2 < sums.size(); c2++) {
            double target_end =
                prec_mass - sums[c2] - water - trunc_shift + prms[start];
            ends.clear();
            findMasses(prms, first_end_bp, len, target_end, prec_tole, &ends);
            for (size_t e = 0; e < ends.size(); e++) {
              addProteoform(start, ends[e] - 1, -1, m, 0.0, trunc_shift);
            }
          }
        }
      }
    }
  }
}

// multi and topmg stages: the protein alone
void TagFilter::geneProtPrsm(const Candidate& cand,
                             const MsHeaderPtr& header_ptr, int group_spec_num,
                             std::vector<SimplePrsmPtrVec>* results) const {
  const FastaSeqPtr& seq_ptr = index_ptr_->getFastaSeqPtr(cand.prot_idx);
  ProteoformPtr proteo_ptr = proteoform_factory::geneDbProteoformPtr(
      seq_ptr, mng_ptr_->prsm_para_ptr_->getFixModPtrVec());
  SimplePrsmPtr prsm_ptr = std::make_shared<SimplePrsm>(
      header_ptr, group_spec_num, proteo_ptr, cand.score);
  (*results)[0].push_back(prsm_ptr);
}

}  // namespace toppic
