/*
 * ps_align.cpp
 *
 *  Created on: Jan 8, 2014
 *      Author: xunlikun
 */
#include "ptmsearch/ps_align.hpp"

namespace prot {

PSAlign::PSAlign() {};

PSAlign::PSAlign(std::vector<double> sp_masses, std::vector<double> seq_masses,
                 BasicDiagPairDiagPtrVec diagonals, PtmMngPtr mng) {
  mng_ = mng;
  sp_masses_ = sp_masses;
  seq_masses_ = seq_masses;
  diagonals_ = diagonals;
}

void PSAlign::compute(SemiAlignTypePtr align_type) {
  initDPPair();
  dp(align_type);
  backtrace();
}

void PSAlign::initDPPair() {
  // init 2d dp pairs
  dp_2d_pairs_.clear();
  segment_bgn_pairs_.clear();
  segment_end_pairs_.clear();
  for (size_t i = 0; i < diagonals_.size(); i++) {
    DPPairPtrVec temp_dppair;
    dp_2d_pairs_.push_back(temp_dppair);
    for (size_t j = 0; j < diagonals_[i]->size(); j++) {
      int x = diagonals_[i]->getDiagPair(j)->getX();
      int y = diagonals_[i]->getDiagPair(j)->getY();
      double score = diagonals_[i]->getDiagPair(j)->getScore();
      double diff = diagonals_[i]->getDiagPair(j)->getDiff();
      dp_2d_pairs_[i].push_back(
          DPPairPtr(
              new DPPair(x, y, score, diff, j, mng_->n_unknown_shift_,
                         diagonals_[i]->getHeader())));
    }
    segment_bgn_pairs_.push_back(dp_2d_pairs_[i][0]);
    segment_end_pairs_.push_back(dp_2d_pairs_[i][diagonals_[i]->size() - 1]);
  }

  // init 1d dp pairs 
  dp_pairs_.clear();
  first_pair_ = DPPairPtr(
      new DPPair(-1, -1, 0, 0, -1, mng_->n_unknown_shift_, nullptr));
  first_pair_->setDiagPrev(nullptr);
  dp_pairs_.push_back(first_pair_);
  for (size_t i = 0; i < dp_2d_pairs_.size(); i++) {
    for (size_t j = 0; j < dp_2d_pairs_[i].size(); j++) {
      dp_pairs_.push_back(dp_2d_pairs_[i][j]);
      if (j > 0) {
        dp_2d_pairs_[i][j]->setDiagPrev(dp_2d_pairs_[i][j - 1]);
      }
    }
  }
  // last pair
  double diff = sp_masses_[sp_masses_.size() - 1]
      - seq_masses_[seq_masses_.size() - 1];
  last_pair_ = DPPairPtr(
      new DPPair(sp_masses_.size(), seq_masses_.size(), 0, diff, -1,
                 mng_->n_unknown_shift_, nullptr));
  last_pair_->setDiagPrev(nullptr);
  dp_pairs_.push_back(last_pair_);
}

void PSAlign::dpPrep() {
  std::sort(dp_pairs_.begin(),dp_pairs_.end(),prot::comparePairUp);
  for (int s = 0; s < mng_->n_unknown_shift_ + 1; s++) {
    dp_pairs_[0]->updateTable(s, 0, PATH_TYPE_NULL, nullptr);
  }
}

DPPairPtr PSAlign::getTruncPre(DPPairPtr cur_pair, int s,
                               SemiAlignTypePtr type) {
  DPPairPtr trunc_prev;
  if (cur_pair == last_pair_) {
    double trunc_score = - std::numeric_limits<double>::max();
    for (size_t i = 0; i < segment_end_pairs_.size(); i++) {
      DPPairPtr prev_pair = segment_end_pairs_[i];
      if (type == SemiAlignTypeFactory::getCompletePtr()
          || type == SemiAlignTypeFactory::getSuffixPtr()) {
        if (prev_pair->getDiagonalHeader()->isProtCTermMatch()) {
          if (prev_pair->getScr(s) > trunc_score) {
            trunc_prev = prev_pair;
            trunc_score = prev_pair->getScr(s);
          }
        }
      }
      else {
        if (prev_pair->getDiagonalHeader()->isPepCTermMatch()) {
          if (prev_pair->getScr(s) > trunc_score) {
            trunc_prev = prev_pair;
            trunc_score = prev_pair->getScr(s);
          }
        }
      }
    }
  } 
  else {
    // if cur_pair is the first in a diagonal 
    if (cur_pair->getDiagOrder() == 0) {
      if (type == SemiAlignTypeFactory::getCompletePtr()
          || type == SemiAlignTypeFactory::getPrefixPtr()) {
        if (cur_pair->getDiagonalHeader()->isProtNTermMatch()) {
          trunc_prev = first_pair_;
        }
      } else {
        if (cur_pair->getDiagonalHeader()->isPepNTermMatch()) {
          trunc_prev = first_pair_;
        }
      }
    }
  }
  return trunc_prev;
}

DPPairPtr PSAlign::getShiftPre(DPPairPtr cur_pair, int p, int s,
                               SemiAlignTypePtr type) {
  if (s < 1) {
    return nullptr;
  }
  // last pair can not shift, only trunc is allowed;
  if (cur_pair == last_pair_) {
    return nullptr;
  }
  int cur_x = cur_pair->getX();
  int cur_y = cur_pair->getY();
  DPPairPtr shift_prev = nullptr;
  double shift_score = -std::numeric_limits<double>::max();

  // cur pair is not last pair
  // first pair can not shift, so q starts from 1
  for (int q = 1; q < p; q++) {
    DPPairPtr prev_pair = dp_pairs_[q];
    int prev_x = prev_pair->getX();
    int prev_y = prev_pair->getY();
    double prev_pair_nterm_shift = prev_pair->getDiagonalHeader()
        ->getProtNTermShift();
    double cur_pair_nterm_shift = cur_pair->getDiagonalHeader()
        ->getProtNTermShift();
    double abs_shift = std::abs(prev_pair_nterm_shift - cur_pair_nterm_shift);
    if (prev_x >= cur_x || prev_y >= cur_y
        || prev_pair->getDiagonalHeader() == cur_pair->getDiagonalHeader()
        || abs_shift > mng_->align_max_shift_) {
      continue;
    }
    double prev_score = prev_pair->getScr(s - 1);
    if (abs_shift > mng_-> align_large_shift_thresh_) {
      prev_score = prev_score - mng_-> align_large_shift_panelty_;
    }
    if (prev_score > shift_score) {
      shift_prev = prev_pair;
      shift_score = prev_score;
    }
  }
  return shift_prev;
}

void PSAlign::dp(SemiAlignTypePtr align_type) {
  dpPrep();
  for (size_t p = 1; p < dp_pairs_.size(); p++) {
    for (int s = 0; s <= mng_->n_unknown_shift_; s++) {
      DPPairPtr trunc_prev = getTruncPre(dp_pairs_[p], s, align_type);
      double trunc_score;
      if (trunc_prev == nullptr) {
        trunc_score = - std::numeric_limits<double>::max();
      } else {
        trunc_score = trunc_prev->getScr(s);
      }

      DPPairPtr diag_prev = dp_pairs_[p]->getDiagPrev();
      double diag_score;
      if (diag_prev != nullptr) {
        diag_score = diag_prev->getScr(s);
      } else {
        diag_score = - std::numeric_limits<double>::max();
      }
      DPPairPtr shift_prev = getShiftPre(dp_pairs_[p], p, s, align_type);
      double shift_score;
      if (shift_prev == nullptr) {
        shift_score =  - std::numeric_limits<double>::max();
      } else {
        shift_score = shift_prev->getScr(s - 1);
      }
      double new_score = dp_pairs_[p]->getPairScore();
      if (trunc_score >= diag_score && trunc_score >= shift_score) {
        if (trunc_score ==  - std::numeric_limits<double>::max()) {
          dp_pairs_[p]->updateTable(s, -std::numeric_limits<double>::max(), 
                                    PATH_TYPE_NULL, nullptr);
        } else {
          dp_pairs_[p]->updateTable(s, trunc_score + new_score, PATH_TYPE_TRUNC,
                                    trunc_prev);
        }
      } else if (diag_score >= shift_score) {
        dp_pairs_[p]->updateTable(s, diag_score + new_score, PATH_TYPE_DIAGONAL,
                                  diag_prev);
      } else {
        dp_pairs_[p]->updateTable(s, shift_score + new_score, PATH_TYPE_SHIFT,
                                  shift_prev);
      }
    }
  }
}

void PSAlign::backtrace() {
  align_scores_.clear();
  backtrack_diagonals_.clear();
  for (int s = 0; s <= mng_->n_unknown_shift_; s++) {
    align_scores_.push_back(last_pair_->getScr(s));
    backtrack_diagonals_.push_back(backtrace(s));
  }
}

DiagonalHeaderPtrVec PSAlign::backtrace(int s) {
  DiagonalHeaderPtrVec list;
  DiagonalHeaderPtr cur_header;
  int cur_end = -1;
  int cur_bgn = -1;
  DPPairPtr p = last_pair_;
  if (p->getPre(s) == nullptr || p->getPre(s) == first_pair_ || p->getScr(s) <= 0) {
    return list;
  }

  while (p != first_pair_) {
    DPPairPtr pre = p->getPre(s);
    if (p == last_pair_) {
      cur_header = pre->getDiagonalHeader();
      cur_end = pre->getY();
    } else if (pre == first_pair_) {
      cur_bgn = p->getY();
      list.push_back(geneDiagonalHeaderPtr(cur_bgn, cur_end, cur_header));
    } else {
      if (p->getType(s) == PATH_TYPE_SHIFT) {
        cur_bgn = p->getY();
        list.push_back(geneDiagonalHeaderPtr(cur_bgn, cur_end, cur_header));
        cur_header = pre->getDiagonalHeader();
        cur_end = pre->getY();
      }
    }
    if (p->getType(s) == PATH_TYPE_SHIFT) {
      s--;
    }
    p = pre;
  }
  std::reverse(list.begin(), list.end());
  return list;
}
} /* namespace prot */
