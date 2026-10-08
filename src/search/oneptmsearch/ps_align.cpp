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

#include "search/oneptmsearch/ps_align.hpp"

#include <algorithm>
#include <limits>
#include <sstream>

#include "search/diag/diag_header_util.hpp"
#include "search/diag/second_align_util.hpp"
#include "search/oneptmsearch/path_type.hpp"

namespace toppic {

PsAlign::PsAlign(const std::vector<double>& ms_masses,
                 const std::vector<double>& seq_masses,
                 const DiagonalPtrVec& diagonal_ptrs,
                 const PsAlignParaPtr& para_ptr)
    : para_ptr_(para_ptr),
      ms_masses_(ms_masses),
      seq_masses_(seq_masses),
      diagonal_ptrs_(diagonal_ptrs) {
  initDpPair();
}

void PsAlign::compute(const ProteoformTypePtr& align_type_ptr) {
  align_type_ptr_ = align_type_ptr;
  dp(align_type_ptr);
  backtrace();
}

void PsAlign::initDpPair() {
  // init 2d dp pairs
  dp_2d_pair_ptrs_.clear();
  segment_bgn_pair_ptrs_.clear();
  segment_end_pair_ptrs_.clear();
  for (size_t i = 0; i < diagonal_ptrs_.size(); i++) {
    DpPairPtrVec temp_dp_pair_ptrs;
    dp_2d_pair_ptrs_.push_back(temp_dp_pair_ptrs);
    for (size_t j = 0; j < diagonal_ptrs_[i]->size(); j++) {
      int x = diagonal_ptrs_[i]->getDiagPair(j)->getX();
      int y = diagonal_ptrs_[i]->getDiagPair(j)->getY();
      double score = diagonal_ptrs_[i]->getDiagPair(j)->getScore();
      double diff = diagonal_ptrs_[i]->getDiagPair(j)->getDiff();
      DpPairPtr dp_pair_ptr = std::make_shared<DpPair>(
          x, y, score, diff, j, para_ptr_->n_unknown_shift_,
          diagonal_ptrs_[i]->getHeader());
      dp_2d_pair_ptrs_[i].push_back(dp_pair_ptr);
    }
    segment_bgn_pair_ptrs_.push_back(dp_2d_pair_ptrs_[i][0]);
    segment_end_pair_ptrs_.push_back(
        dp_2d_pair_ptrs_[i][diagonal_ptrs_[i]->size() - 1]);
  }

  // init 1d dp pairs
  dp_pair_ptrs_.clear();
  first_pair_ptr_ = std::make_shared<DpPair>(
      -1, -1, 0, 0, -1, para_ptr_->n_unknown_shift_, nullptr);
  first_pair_ptr_->setDiagPrevPairPtr(nullptr);
  dp_pair_ptrs_.push_back(first_pair_ptr_);
  for (size_t i = 0; i < dp_2d_pair_ptrs_.size(); i++) {
    for (size_t j = 0; j < dp_2d_pair_ptrs_[i].size(); j++) {
      dp_pair_ptrs_.push_back(dp_2d_pair_ptrs_[i][j]);
      if (j > 0) {
        dp_2d_pair_ptrs_[i][j]->setDiagPrevPairPtr(dp_2d_pair_ptrs_[i][j - 1]);
      }
    }
  }
  // last pair
  double diff =
      ms_masses_[ms_masses_.size() - 1] - seq_masses_[seq_masses_.size() - 1];
  last_pair_ptr_ =
      std::make_shared<DpPair>(ms_masses_.size(), seq_masses_.size(), 0, diff,
                               -1, para_ptr_->n_unknown_shift_, nullptr);
  last_pair_ptr_->setDiagPrevPairPtr(nullptr);
  dp_pair_ptrs_.push_back(last_pair_ptr_);
  std::sort(dp_pair_ptrs_.begin(), dp_pair_ptrs_.end(), Pair::cmpPosInc);
  // init indexes. The scan over all previous pairs for every pair is the
  // hot loop of the multiple-shift search; it reads only the position and
  // the diagonal of each pair and a per-diagonal-pair shift test, so those
  // are flattened into plain arrays first. The shift between two pairs
  // depends only on their diagonals, not on the pairs, and calling the
  // shared_ptr getters inside the loop cost four atomic reference count
  // updates per visited pair.
  size_t pair_num = dp_pair_ptrs_.size();
  size_t diag_num = diagonal_ptrs_.size();
  std::vector<int> pair_xs(pair_num);
  std::vector<int> pair_ys(pair_num);
  // the diagonal of each pair; -1 for the first and last pair
  std::vector<int> pair_diags(pair_num);
  for (size_t i = 0; i < pair_num; i++) {
    const DpPairPtr& pair_ptr = dp_pair_ptrs_[i];
    pair_xs[i] = pair_ptr->getX();
    pair_ys[i] = pair_ptr->getY();
    if (pair_ptr == first_pair_ptr_ || pair_ptr == last_pair_ptr_) {
      pair_diags[i] = -1;
    } else {
      pair_diags[i] = pair_ptr->getDiagHeader()->getId();
    }
  }
  // shift_allowed[prev][cur]: the shift from diagonal prev to diagonal cur
  // is inside [align_min_shift_, align_max_shift_] and the diagonals differ;
  // shift_large[prev][cur]: its absolute value is above the penalty
  // threshold
  std::vector<double> diag_shifts(diag_num);
  for (size_t d = 0; d < diag_num; d++) {
    diag_shifts[d] = diagonal_ptrs_[d]->getHeader()->getProtNTermShift();
  }
  std::vector<std::vector<char>> shift_allowed(
      diag_num, std::vector<char>(diag_num, 0));
  std::vector<std::vector<char>> shift_large(
      diag_num, std::vector<char>(diag_num, 0));
  for (size_t prev = 0; prev < diag_num; prev++) {
    for (size_t cur = 0; cur < diag_num; cur++) {
      if (prev == cur) {
        continue;
      }
      double shift = diag_shifts[cur] - diag_shifts[prev];
      if (shift > para_ptr_->align_max_shift_ ||
          shift < para_ptr_->align_min_shift_) {
        continue;
      }
      shift_allowed[prev][cur] = 1;
      if (std::abs(shift) > para_ptr_->align_large_shift_thresh_) {
        shift_large[prev][cur] = 1;
      }
    }
  }
  for (size_t i = 0; i < pair_num; i++) {
    std::vector<int> ends(diag_num, -1);
    std::vector<bool> large_shifts(diag_num, false);
    int cur_diag = pair_diags[i];
    if (cur_diag < 0) {
      // the first and last pair
      idxes_.push_back(ends);
      penalties_.push_back(large_shifts);
      continue;
    }
    int cur_x = pair_xs[i];
    int cur_y = pair_ys[i];
    // cur pair is not last pair
    // first pair can not shift, so j starts from 1
    for (size_t j = 1; j < i; j++) {
      int prev_diag = pair_diags[j];
      if (prev_diag < 0 || pair_xs[j] >= cur_x || pair_ys[j] >= cur_y ||
          !shift_allowed[prev_diag][cur_diag]) {
        continue;
      }
      ends[prev_diag] = j;
      if (shift_large[prev_diag][cur_diag]) {
        large_shifts[prev_diag] = true;
      }
    }
    idxes_.push_back(ends);
    penalties_.push_back(large_shifts);
  }
}

void PsAlign::dpPrep() {
  std::sort(dp_pair_ptrs_.begin(), dp_pair_ptrs_.end(), Pair::cmpPosInc);
  for (int s = 0; s < para_ptr_->n_unknown_shift_ + 1; s++) {
    dp_pair_ptrs_[0]->updateTable(s, 0, path_type::TYPE_NULL, nullptr);
  }
}

inline DpPairPtr PsAlign::getTruncPre(const DpPairPtr& cur_pair_ptr, int s,
                                      const ProteoformTypePtr& align_type_ptr) {
  DpPairPtr trunc_prev_ptr;
  if (cur_pair_ptr == last_pair_ptr_) {
    double trunc_score = -std::numeric_limits<double>::max();
    for (size_t i = 0; i < segment_end_pair_ptrs_.size(); i++) {
      DpPairPtr prev_pair_ptr = segment_end_pair_ptrs_[i];
      if (align_type_ptr == ProteoformType::COMPLETE ||
          align_type_ptr == ProteoformType::SUFFIX) {
        if (prev_pair_ptr->getDiagHeader()->isProtCTermMatch()) {
          if (prev_pair_ptr->getScore(s) > trunc_score) {
            trunc_prev_ptr = prev_pair_ptr;
            trunc_score = prev_pair_ptr->getScore(s);
          }
        }
      } else {
        if (prev_pair_ptr->getDiagHeader()->isPepCTermMatch()) {
          if (prev_pair_ptr->getScore(s) > trunc_score) {
            trunc_prev_ptr = prev_pair_ptr;
            trunc_score = prev_pair_ptr->getScore(s);
          }
        }
      }
    }
  } else {
    // if cur_pair_ptr is the first in a diagonal
    if (cur_pair_ptr->getDiagOrder() == 0) {
      if (align_type_ptr == ProteoformType::COMPLETE ||
          align_type_ptr == ProteoformType::PREFIX) {
        if (cur_pair_ptr->getDiagHeader()->isProtNTermMatch()) {
          trunc_prev_ptr = first_pair_ptr_;
        }
      } else {
        if (cur_pair_ptr->getDiagHeader()->isPepNTermMatch()) {
          trunc_prev_ptr = first_pair_ptr_;
        }
      }
    }
  }
  return trunc_prev_ptr;
}

DpPairPtr PsAlign::getShiftPre(int p, int s, double& shift_score) {
  shift_score = -std::numeric_limits<double>::max();
  if (s < 1) {
    return nullptr;
  }
  // last pair can not shift, only trunc is allowed;
  DpPairPtr shift_prev = nullptr;

  // cur pair is not last pair
  // first pair can not shift, so q starts from 1
  for (size_t d = 0; d < idxes_[p].size(); d++) {
    if (idxes_[p][d] >= 0) {
      DpPairPtr prev_pair_ptr = dp_pair_ptrs_[idxes_[p][d]];
      double prev_score = prev_pair_ptr->getScore(s - 1);
      if (penalties_[p][d]) {
        prev_score = prev_score - para_ptr_->align_large_shift_panelty_;
      }
      if (prev_score > shift_score) {
        shift_prev = prev_pair_ptr;
        shift_score = prev_score;
      }
    }
  }
  return shift_prev;
}

void PsAlign::dp(const ProteoformTypePtr& align_type_ptr) {
  dpPrep();
  for (size_t p = 1; p < dp_pair_ptrs_.size(); p++) {
    for (int s = 0; s <= para_ptr_->n_unknown_shift_; s++) {
      DpPairPtr trunc_prev_ptr =
          getTruncPre(dp_pair_ptrs_[p], s, align_type_ptr);
      double trunc_score;
      if (trunc_prev_ptr == nullptr) {
        trunc_score = -std::numeric_limits<double>::max();
      } else {
        trunc_score = trunc_prev_ptr->getScore(s);
      }

      DpPairPtr diag_prev = dp_pair_ptrs_[p]->getDiagPrevPairPtr();
      double diag_score;
      if (diag_prev != nullptr) {
        diag_score = diag_prev->getScore(s);
      } else {
        diag_score = -std::numeric_limits<double>::max();
      }
      // The shift predecessor's score with the large-shift penalty taken
      // off, as PsAlignMultiType stores it. The penalty used to decide only
      // which diagonal to shift from, and the stored score was the
      // predecessor's without it (toppic_test experiment 145).
      double shift_score;
      DpPairPtr shift_prev = getShiftPre(p, s, shift_score);
      double new_score = dp_pair_ptrs_[p]->getPairScore();
      if (trunc_score >= diag_score && trunc_score >= shift_score) {
        if (trunc_score == -std::numeric_limits<double>::max()) {
          dp_pair_ptrs_[p]->updateTable(s, -std::numeric_limits<double>::max(),
                                        path_type::TYPE_NULL, nullptr);
        } else {
          dp_pair_ptrs_[p]->updateTable(s, trunc_score + new_score,
                                        path_type::TYPE_TRUNC, trunc_prev_ptr);
        }
      } else if (diag_score >= shift_score) {
        dp_pair_ptrs_[p]->updateTable(s, diag_score + new_score,
                                      path_type::TYPE_DIAGONAL, diag_prev);
      } else {
        dp_pair_ptrs_[p]->updateTable(s, shift_score + new_score,
                                      path_type::TYPE_SHIFT, shift_prev);
      }
    }
  }
}

void PsAlign::backtrace() {
  align_scores_.clear();
  backtrack_diagonal_ptrs_.clear();
  for (int s = 0; s <= para_ptr_->n_unknown_shift_; s++) {
    align_scores_.push_back(last_pair_ptr_->getScore(s));
    backtrack_diagonal_ptrs_.push_back(backtrace(s));
  }
}

DiagHeaderPtrVec PsAlign::backtrace(int s) {
  DiagHeaderPtrVec list;
  DiagHeaderPtr cur_header;
  int cur_end = -1;
  int cur_bgn = -1;
  DpPairPtr p = last_pair_ptr_;
  if (p->getPrevPairPtr(s) == nullptr ||
      p->getPrevPairPtr(s) == first_pair_ptr_ || p->getScore(s) <= 0) {
    return list;
  }

  // LOG_DEBUG("backtrace start ");
  while (p != first_pair_ptr_) {
    DpPairPtr pre = p->getPrevPairPtr(s);
    // LOG_DEBUG("p " << p << " pre ptr " << pre);
    if (p == last_pair_ptr_) {
      cur_header = pre->getDiagHeader();
      cur_end = pre->getY();
    } else if (pre == first_pair_ptr_) {
      cur_bgn = p->getY();
      list.push_back(
          diag_header_util::geneDiagHeaderPtr(cur_bgn, cur_end, cur_header));
    } else {
      if (p->getType(s) == path_type::TYPE_SHIFT) {
        cur_bgn = p->getY();
        list.push_back(
            diag_header_util::geneDiagHeaderPtr(cur_bgn, cur_end, cur_header));
        cur_header = pre->getDiagHeader();
        cur_end = pre->getY();
      }
    }
    if (p->getType(s) == path_type::TYPE_SHIFT) {
      s--;
    }
    p = pre;
  }
  std::reverse(list.begin(), list.end());
  return list;
}

PrsmPtr PsAlign::geneResult(int shift_num, const ProteoformPtr& proteo_ptr,
                            DeconvMsPtrVec& deconv_ms_ptr_vec,
                            ExtendMsPtrVec& ms_three_ptr_vec,
                            const PrsmParaPtr& prsm_para_ptr) {
  DiagHeaderPtrVec header_ptrs = getDiagHeaders(shift_num);
  if (header_ptrs.size() == 0) {
    return nullptr;
  }
  // the first alignment's PrSM: main's result
  double refine_prec_mass = 0;
  ExtendMsPtrVec refine_ms_ptr_vec;
  PrsmPtr first_prsm_ptr = second_align_util::genePrsm(
      proteo_ptr, header_ptrs, deconv_ms_ptr_vec, ms_three_ptr_vec,
      prsm_para_ptr, refine_prec_mass, refine_ms_ptr_vec);
  if (!prsm_para_ptr->useSecondAlign() || first_prsm_ptr == nullptr) {
    return first_prsm_ptr;
  }

  // the second alignment at the refined precursor
  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  double prec_mass =
      ms_three_ptr_vec[0]->getMsHeaderPtr()->getFirstPrecMonoMass();
  PrmPeakPtrVec prm_peaks = second_align_util::genePrmPeaks(
      deconv_ms_ptr_vec, sp_para_ptr, refine_prec_mass);
  DiagHeaderPtrVec moved_header_ptrs = second_align_util::moveHeaders(
      diagonal_ptrs_, refine_prec_mass - prec_mass, refine_prec_mass,
      proteo_ptr);
  std::vector<double> ms_masses;
  DiagonalPtrVec diagonal_ptrs = second_align_util::geneDiagonals(
      moved_header_ptrs, prm_peaks, deconv_ms_ptr_vec.size(), proteo_ptr,
      ms_masses);
  PsAlign second_align(ms_masses, seq_masses_, diagonal_ptrs, para_ptr_);
  second_align.compute(align_type_ptr_);
  DiagHeaderPtrVec second_header_ptrs = second_align.getDiagHeaders(shift_num);

  PrsmPtr second_prsm_ptr = nullptr;
  double second_prec_mass = 0;
  if (second_header_ptrs.size() > 0) {
    ExtendMsPtrVec second_ms_ptr_vec;
    second_prsm_ptr = second_align_util::genePrsm(
        proteo_ptr, second_header_ptrs, deconv_ms_ptr_vec, refine_ms_ptr_vec,
        prsm_para_ptr, second_prec_mass, second_ms_ptr_vec, true);
  }
  if (second_align_util::logEnabled()) {
    std::ostringstream out;
    out.precision(4);
    out << std::fixed;
    out << "one\t" << deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getSpecId()
        << "\t" << align_type_ptr_->getName() << "\t" << shift_num << "\t"
        << proteo_ptr->getSeqName() << "\t" << prec_mass << "\t"
        << getAlignScr(shift_num) << "\t" << refine_prec_mass << "\t"
        << second_align_util::prsmInfo(first_prsm_ptr) << "\t";
    if (second_header_ptrs.size() > 0) {
      out << second_align.getAlignScr(shift_num) << "\t" << second_prec_mass;
    } else {
      out << "-\t-";
    }
    out << "\t" << second_align_util::prsmInfo(second_prsm_ptr) << "\t"
        << (second_prsm_ptr == nullptr ? "first" : "second");
    second_align_util::log(out.str());
  }
  if (second_prsm_ptr == nullptr) {
    return first_prsm_ptr;
  }
  return second_prsm_ptr;
}

} /* namespace toppic */
