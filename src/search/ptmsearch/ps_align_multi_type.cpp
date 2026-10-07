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

#include "search/ptmsearch/ps_align_multi_type.hpp"

#include <algorithm>
#include <limits>
#include <sstream>

#include "search/diag/diag_header_util.hpp"
#include "search/diag/second_align_util.hpp"
#include "search/oneptmsearch/path_type.hpp"

namespace toppic {

namespace {

// The four proteoform types differ only in which terminus the alignment is
// required to start/end at (search/oneptmsearch/ps_align.cpp getTruncPre):
// COMPLETE and PREFIX must start at the protein N-terminus, SUFFIX and
// INTERNAL may start anywhere; COMPLETE and SUFFIX must end at the protein
// C-terminus, PREFIX and INTERNAL may end anywhere. These two predicates are
// the only type-dependent part of the dynamic program.
inline bool startUsesProtNTerm(int type_id) {
  return type_id == ProteoformType::COMPLETE->getId() ||
         type_id == ProteoformType::PREFIX->getId();
}

inline bool endUsesProtCTerm(int type_id) {
  return type_id == ProteoformType::COMPLETE->getId() ||
         type_id == ProteoformType::SUFFIX->getId();
}

}  // namespace

PsAlignMultiType::PsAlignMultiType(const std::vector<double>& ms_masses,
                                   const std::vector<double>& seq_masses,
                                   const DiagonalPtrVec& diagonal_ptrs,
                                   const PsAlignParaPtr& para_ptr)
    : para_ptr_(para_ptr),
      ms_masses_(ms_masses),
      seq_masses_(seq_masses),
      diagonal_ptrs_(diagonal_ptrs) {
  align_scores_.resize(kProteoformTypeNum);
  backtrack_diagonal_ptrs_.resize(kProteoformTypeNum);
  initDpPair();
}

void PsAlignMultiType::compute(const ProteoformTypePtrVec& type_ptrs) {
  type_ptrs_ = type_ptrs;
  second_align_ptrs_.clear();
  dp(type_ptrs);
  backtrace(type_ptrs);
}

void PsAlignMultiType::initDpPair() {
  // init 2d dp pairs
  dp_2d_pair_ptrs_.clear();
  segment_bgn_pair_ptrs_.clear();
  segment_end_pair_ptrs_.clear();
  for (size_t i = 0; i < diagonal_ptrs_.size(); i++) {
    MultiTypeDpPairPtrVec temp_dp_pair_ptrs;
    dp_2d_pair_ptrs_.push_back(temp_dp_pair_ptrs);
    for (size_t j = 0; j < diagonal_ptrs_[i]->size(); j++) {
      int x = diagonal_ptrs_[i]->getDiagPair(j)->getX();
      int y = diagonal_ptrs_[i]->getDiagPair(j)->getY();
      double score = diagonal_ptrs_[i]->getDiagPair(j)->getScore();
      double diff = diagonal_ptrs_[i]->getDiagPair(j)->getDiff();
      MultiTypeDpPairPtr dp_pair_ptr = std::make_shared<MultiTypeDpPair>(
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
  first_pair_ptr_ = std::make_shared<MultiTypeDpPair>(
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
  last_pair_ptr_ = std::make_shared<MultiTypeDpPair>(
      ms_masses_.size(), seq_masses_.size(), 0, diff, -1,
      para_ptr_->n_unknown_shift_, nullptr);
  last_pair_ptr_->setDiagPrevPairPtr(nullptr);
  dp_pair_ptrs_.push_back(last_pair_ptr_);
  std::sort(dp_pair_ptrs_.begin(), dp_pair_ptrs_.end(), Pair::cmpPosInc);
  // init indexes. Same precomputation as PsAlign::initDpPair: it is
  // type-independent (it depends only on pair position and diagonal), so it
  // is still done exactly once here and shared by every type computed by
  // this object, the same as it is already shared across PsAlign::compute's
  // repeated calls for different types.
  size_t pair_num = dp_pair_ptrs_.size();
  size_t diag_num = diagonal_ptrs_.size();
  std::vector<int> pair_xs(pair_num);
  std::vector<int> pair_ys(pair_num);
  std::vector<int> pair_diags(pair_num);
  for (size_t i = 0; i < pair_num; i++) {
    const MultiTypeDpPairPtr& pair_ptr = dp_pair_ptrs_[i];
    pair_xs[i] = pair_ptr->getX();
    pair_ys[i] = pair_ptr->getY();
    if (pair_ptr == first_pair_ptr_ || pair_ptr == last_pair_ptr_) {
      pair_diags[i] = -1;
    } else {
      pair_diags[i] = pair_ptr->getDiagHeader()->getId();
    }
  }
  std::vector<double> diag_shifts(diag_num);
  for (size_t d = 0; d < diag_num; d++) {
    diag_shifts[d] = diagonal_ptrs_[d]->getHeader()->getProtNTermShift();
  }
  std::vector<std::vector<char>> shift_allowed(diag_num,
                                               std::vector<char>(diag_num, 0));
  std::vector<std::vector<char>> shift_large(diag_num,
                                             std::vector<char>(diag_num, 0));
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
      idxes_.push_back(ends);
      penalties_.push_back(large_shifts);
      continue;
    }
    int cur_x = pair_xs[i];
    int cur_y = pair_ys[i];
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

void PsAlignMultiType::dpPrep(const ProteoformTypePtrVec& type_ptrs) {
  // dp_pair_ptrs_ is sorted once, in initDpPair(); the position order of the
  // pairs never changes between compute() calls, so unlike PsAlign::dpPrep
  // it is not re-sorted here.
  for (size_t ti = 0; ti < type_ptrs.size(); ti++) {
    int t = type_ptrs[ti]->getId();
    for (int s = 0; s <= para_ptr_->n_unknown_shift_; s++) {
      dp_pair_ptrs_[0]->updateTable(t, s, 0, path_type::TYPE_NULL, nullptr);
    }
  }
}

void PsAlignMultiType::getTruncPreLast(
    const ProteoformTypePtrVec& type_ptrs, int s,
    std::array<MultiTypeDpPairPtr, kProteoformTypeNum>* trunc_prev,
    std::array<double, kProteoformTypeNum>* trunc_score) {
  double neg_inf = -std::numeric_limits<double>::max();
  for (size_t ti = 0; ti < type_ptrs.size(); ti++) {
    int t = type_ptrs[ti]->getId();
    (*trunc_prev)[t] = nullptr;
    (*trunc_score)[t] = neg_inf;
  }
  for (size_t i = 0; i < segment_end_pair_ptrs_.size(); i++) {
    const MultiTypeDpPairPtr& prev_pair_ptr = segment_end_pair_ptrs_[i];
    DiagHeaderPtr header_ptr = prev_pair_ptr->getDiagHeader();
    bool prot_c_term = header_ptr->isProtCTermMatch();
    bool pep_c_term = header_ptr->isPepCTermMatch();
    for (size_t ti = 0; ti < type_ptrs.size(); ti++) {
      int t = type_ptrs[ti]->getId();
      bool ok = endUsesProtCTerm(t) ? prot_c_term : pep_c_term;
      if (!ok) {
        continue;
      }
      double score = prev_pair_ptr->getScore(t, s);
      if (score > (*trunc_score)[t]) {
        (*trunc_score)[t] = score;
        (*trunc_prev)[t] = prev_pair_ptr;
      }
    }
  }
}

void PsAlignMultiType::dp(const ProteoformTypePtrVec& type_ptrs) {
  dpPrep(type_ptrs);
  double neg_inf = -std::numeric_limits<double>::max();
  size_t last_p = dp_pair_ptrs_.size() - 1;
  for (size_t p = 1; p < dp_pair_ptrs_.size(); p++) {
    const MultiTypeDpPairPtr& cur_pair_ptr = dp_pair_ptrs_[p];
    MultiTypeDpPairPtr diag_prev = cur_pair_ptr->getDiagPrevPairPtr();
    double new_score = cur_pair_ptr->getPairScore();
    bool is_last = (p == last_p);
    // Whether cur_pair_ptr is the first pair of its diagonal segment (the
    // only case, besides the last virtual pair, where a truncation
    // transition is possible); the two start-truncation flags do not depend
    // on type, so they are read once per p and reused by every type below,
    // exactly as PsAlign::getTruncPre reads them once per (p, type) call.
    bool is_diag_start = (!is_last) && (cur_pair_ptr->getDiagOrder() == 0);
    bool start_prot_ok = false;
    bool start_pep_ok = false;
    if (is_diag_start) {
      start_prot_ok = cur_pair_ptr->getDiagHeader()->isProtNTermMatch();
      start_pep_ok = cur_pair_ptr->getDiagHeader()->isPepNTermMatch();
    }

    const std::vector<int>& ends = idxes_[p];
    const std::vector<bool>& large = penalties_[p];

    for (int s = 0; s <= para_ptr_->n_unknown_shift_; s++) {
      // Shift transition: the candidate predecessors (one per diagonal) are
      // the same for every type, so the scan over diagonals is done once
      // per (p, s) and updates every requested type's running best shift
      // predecessor as it goes, instead of PsAlign::getShiftPre's one full
      // scan per (p, s, type).
      std::array<MultiTypeDpPairPtr, kProteoformTypeNum> shift_prev{};
      std::array<double, kProteoformTypeNum> shift_score;
      shift_score.fill(neg_inf);
      if (s >= 1) {
        for (size_t d = 0; d < ends.size(); d++) {
          int idx = ends[d];
          if (idx < 0) {
            continue;
          }
          const MultiTypeDpPairPtr& prev_pair_ptr = dp_pair_ptrs_[idx];
          bool penalty = large[d];
          for (size_t ti = 0; ti < type_ptrs.size(); ti++) {
            int t = type_ptrs[ti]->getId();
            double prev_score = prev_pair_ptr->getScore(t, s - 1);
            if (penalty) {
              prev_score -= para_ptr_->align_large_shift_panelty_;
            }
            if (prev_score > shift_score[t]) {
              shift_score[t] = prev_score;
              shift_prev[t] = prev_pair_ptr;
            }
          }
        }
      }

      // Truncation transition at the last (virtual) pair: likewise fused
      // across types into one scan of segment_end_pair_ptrs_.
      std::array<MultiTypeDpPairPtr, kProteoformTypeNum> trunc_prev{};
      std::array<double, kProteoformTypeNum> trunc_score;
      trunc_score.fill(neg_inf);
      if (is_last) {
        getTruncPreLast(type_ptrs, s, &trunc_prev, &trunc_score);
      }

      for (size_t ti = 0; ti < type_ptrs.size(); ti++) {
        const ProteoformTypePtr& type_ptr = type_ptrs[ti];
        int t = type_ptr->getId();

        MultiTypeDpPairPtr cur_trunc_prev;
        double cur_trunc_score = neg_inf;
        if (is_last) {
          cur_trunc_prev = trunc_prev[t];
          cur_trunc_score = trunc_score[t];
        } else if (is_diag_start) {
          bool ok = startUsesProtNTerm(t) ? start_prot_ok : start_pep_ok;
          if (ok) {
            cur_trunc_prev = first_pair_ptr_;
          }
        }
        if (cur_trunc_prev == nullptr) {
          cur_trunc_score = neg_inf;
        } else if (!is_last) {
          cur_trunc_score = cur_trunc_prev->getScore(t, s);
        }

        double diag_score =
            (diag_prev != nullptr) ? diag_prev->getScore(t, s) : neg_inf;
        double cur_shift_score = shift_score[t];
        MultiTypeDpPairPtr cur_shift_prev = shift_prev[t];

        if (cur_trunc_score >= diag_score &&
            cur_trunc_score >= cur_shift_score) {
          if (cur_trunc_score == neg_inf) {
            cur_pair_ptr->updateTable(t, s, neg_inf, path_type::TYPE_NULL,
                                      nullptr);
          } else {
            cur_pair_ptr->updateTable(t, s, cur_trunc_score + new_score,
                                      path_type::TYPE_TRUNC, cur_trunc_prev);
          }
        } else if (diag_score >= cur_shift_score) {
          cur_pair_ptr->updateTable(t, s, diag_score + new_score,
                                    path_type::TYPE_DIAGONAL, diag_prev);
        } else {
          cur_pair_ptr->updateTable(t, s, cur_shift_score + new_score,
                                    path_type::TYPE_SHIFT, cur_shift_prev);
        }
      }
    }
  }
}

void PsAlignMultiType::backtrace(const ProteoformTypePtrVec& type_ptrs) {
  for (size_t ti = 0; ti < type_ptrs.size(); ti++) {
    int t = type_ptrs[ti]->getId();
    align_scores_[t].clear();
    backtrack_diagonal_ptrs_[t].clear();
    for (int s = 0; s <= para_ptr_->n_unknown_shift_; s++) {
      align_scores_[t].push_back(last_pair_ptr_->getScore(t, s));
      backtrack_diagonal_ptrs_[t].push_back(backtrace(t, s));
    }
  }
}

DiagHeaderPtrVec PsAlignMultiType::backtrace(int type_id, int s) {
  DiagHeaderPtrVec list;
  DiagHeaderPtr cur_header;
  int cur_end = -1;
  int cur_bgn = -1;
  MultiTypeDpPairPtr p = last_pair_ptr_;
  if (p->getPrevPairPtr(type_id, s) == nullptr ||
      p->getPrevPairPtr(type_id, s) == first_pair_ptr_ ||
      p->getScore(type_id, s) <= 0) {
    return list;
  }

  while (p != first_pair_ptr_) {
    MultiTypeDpPairPtr pre = p->getPrevPairPtr(type_id, s);
    if (p == last_pair_ptr_) {
      cur_header = pre->getDiagHeader();
      cur_end = pre->getY();
    } else if (pre == first_pair_ptr_) {
      cur_bgn = p->getY();
      list.push_back(
          diag_header_util::geneDiagHeaderPtr(cur_bgn, cur_end, cur_header));
    } else {
      if (p->getType(type_id, s) == path_type::TYPE_SHIFT) {
        cur_bgn = p->getY();
        list.push_back(
            diag_header_util::geneDiagHeaderPtr(cur_bgn, cur_end, cur_header));
        cur_header = pre->getDiagHeader();
        cur_end = pre->getY();
      }
    }
    if (p->getType(type_id, s) == path_type::TYPE_SHIFT) {
      s--;
    }
    p = pre;
  }
  std::reverse(list.begin(), list.end());
  return list;
}

PrsmPtr PsAlignMultiType::geneResult(const ProteoformTypePtr& type_ptr,
                                     int shift_num,
                                     const ProteoformPtr& proteo_ptr,
                                     DeconvMsPtrVec& deconv_ms_ptr_vec,
                                     ExtendMsPtrVec& ms_three_ptr_vec,
                                     const PrsmParaPtr& prsm_para_ptr) {
  DiagHeaderPtrVec header_ptrs = getDiagHeaders(type_ptr, shift_num);
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

  // the second alignment at the refined precursor, shared by the types
  // whose first results refine the precursor to the same mass
  double prec_mass =
      ms_three_ptr_vec[0]->getMsHeaderPtr()->getFirstPrecMonoMass();
  std::shared_ptr<PsAlignMultiType> second_align_ptr;
  for (size_t i = 0; i < second_align_ptrs_.size(); i++) {
    if (second_align_ptrs_[i].first == refine_prec_mass) {
      second_align_ptr = second_align_ptrs_[i].second;
      break;
    }
  }
  if (second_align_ptr == nullptr) {
    SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
    bool strict = prsm_para_ptr->getSecondAlignMode() == "strict";
    PrmPeakPtrVec prm_peaks = second_align_util::genePrmPeaks(
        deconv_ms_ptr_vec, sp_para_ptr, refine_prec_mass, strict);
    DiagHeaderPtrVec moved_header_ptrs = second_align_util::moveHeaders(
        diagonal_ptrs_, refine_prec_mass - prec_mass, refine_prec_mass,
        proteo_ptr);
    std::vector<double> ms_masses;
    DiagonalPtrVec diagonal_ptrs = second_align_util::geneDiagonals(
        moved_header_ptrs, prm_peaks, deconv_ms_ptr_vec.size(), proteo_ptr,
        ms_masses);
    second_align_ptr = std::make_shared<PsAlignMultiType>(
        ms_masses, seq_masses_, diagonal_ptrs, para_ptr_);
    second_align_ptr->compute(type_ptrs_);
    second_align_ptrs_.push_back(
        std::make_pair(refine_prec_mass, second_align_ptr));
  }
  DiagHeaderPtrVec second_header_ptrs =
      second_align_ptr->getDiagHeaders(type_ptr, shift_num);

  PrsmPtr second_prsm_ptr = nullptr;
  double second_prec_mass = 0;
  if (second_header_ptrs.size() > 0) {
    ExtendMsPtrVec second_ms_ptr_vec;
    second_prsm_ptr = second_align_util::genePrsm(
        proteo_ptr, second_header_ptrs, deconv_ms_ptr_vec, refine_ms_ptr_vec,
        prsm_para_ptr, second_prec_mass, second_ms_ptr_vec,
        prsm_para_ptr->useSecondAlignClamp() ? prec_mass : 0);
  }
  if (second_align_util::logEnabled()) {
    std::ostringstream out;
    out.precision(4);
    out << std::fixed;
    out << "multi\t" << deconv_ms_ptr_vec[0]->getMsHeaderPtr()->getSpecId()
        << "\t" << type_ptr->getName() << "\t" << shift_num << "\t"
        << proteo_ptr->getSeqName() << "\t" << prec_mass << "\t"
        << getAlignScr(type_ptr, shift_num) << "\t"
        << refine_prec_mass << "\t"
        << second_align_util::prsmInfo(first_prsm_ptr) << "\t";
    if (second_header_ptrs.size() > 0) {
      out << second_align_ptr->getAlignScr(type_ptr, shift_num) << "\t"
          << second_prec_mass;
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
