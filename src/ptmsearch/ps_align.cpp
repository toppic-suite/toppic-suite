/*
 * ps_align.cpp
 *
 *  Created on: Jan 8, 2014
 *      Author: xunlikun
 */
#include "base/semi_align_type.hpp"
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
  for (unsigned int i = 0; i < diagonals_.size(); i++) {
    DPPairPtrVec temp_dppair;
    dp_2d_pairs_.push_back(temp_dppair);
    for (unsigned int j = 0; j < diagonals_[i]->size(); j++) {
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
  for (unsigned int i = 0; i < dp_2d_pairs_.size(); i++) {
    for (unsigned int j = 0; j < dp_2d_pairs_[i].size(); j++) {
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
    for (unsigned int i = 0; i < segment_end_pairs_.size(); i++) {
      DPPairPtr prev_pair = segment_end_pairs_[i];
      if (type == SemiAlignTypeFactory::getCompletePtr()
          || type == SemiAlignTypeFactory::getSuffixPtr()) {
        if (prev_pair->getDiagonalHeader()->isProtCTermMatch()) {
          if (prev_pair->getSrc(s) > trunc_score) {
            trunc_prev = prev_pair;
            trunc_score = prev_pair->getSrc(s);
          }
        }
      }
      else {
        if (prev_pair->getDiagonalHeader()->isPepCTermMatch()) {
          if (prev_pair->getSrc(s) > trunc_score) {
            trunc_prev = prev_pair;
            trunc_score = prev_pair->getSrc(s);
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
  int cur_x = cur_pair->getX();
  int cur_y = cur_pair->getY();
  DPPairPtr shift_prev = nullptr;
  double shift_score = -std::numeric_limits<double>::max();
  // last pair
  if (cur_pair == last_pair_) {
    double prev_pair_min_n_term_shift = last_pair_->getDiff() - mng_->align_max_shift_;
    double prev_pair_max_n_term_shift = sp_masses_[sp_masses_.size() - 1];
    for (int q = 0; q < p; q++) {
      DPPairPtr prev_pair = dp_pairs_[q];
      int prev_x = prev_pair->getX();
      int prev_y = prev_pair->getY();
      if (prev_x >= cur_x || prev_y >= cur_y
          || prev_pair->getDiagonalHeader()
          == cur_pair->getDiagonalHeader()) {
        continue;
      }
      if ((type == SemiAlignTypeFactory::getCompletePtr()
           || type == SemiAlignTypeFactory::getSuffixPtr())
          && !prev_pair->getDiagonalHeader()->isAlignSuffix()) {
        continue;
      }
      if (type == SemiAlignTypeFactory::getPrefixPtr()
           || type == SemiAlignTypeFactory::getInternalPtr()) {
        double prev_pair_n_term_shift = prev_pair->getDiagonalHeader()
          ->getProtNTermShift();
        if (prev_pair_n_term_shift < prev_pair_min_n_term_shift ||
            prev_pair_n_term_shift > prev_pair_max_n_term_shift) {
          continue;
        }
        double abs_shift = std::abs(prev_pair_n_term_shift - last_pair_->getDiff()); 
        if (abs_shift < mng_->align_min_shift_) {
          continue;
        }
      }
      /*
      if (s == 2 && type == SemiAlignTypeFactory::getSuffixPtr()) {
        std::cout << "prev pair " << " align suffix " << prev_pair->getDiagonalHeader()->isAlignSuffix()
            << " prev n term shift " << prev_pair->getDiagonalHeader()->getProtNTermShift()
            << " last shift " << last_pair_->getDiff()
            << " score " << prev_pair->getSrc(s-1) 
            << std::endl;
      }
      */

      if (prev_pair->getSrc(s - 1) > shift_score) {
        shift_prev = prev_pair;
        shift_score = dp_pairs_[q]->getSrc(s - 1);
      }
    }
  } 
  // not last pair
  else {
    // if complete or prefix alignment, consider first pair only if 
    // isAlignPrefix is true
    if (type == SemiAlignTypeFactory::getCompletePtr()
        || type == SemiAlignTypeFactory::getPrefixPtr()) {
      if (first_pair_->getSrc(s - 1) > shift_score
          && cur_pair->getDiagonalHeader()->isAlignPrefix()) {
        shift_prev = first_pair_;
        shift_score = first_pair_->getSrc(s - 1);
      }
    } 
    else {
      double cur_pair_min_n_term_shift = - seq_masses_[seq_masses_.size()-1];
      double cur_pair_max_n_term_shift = mng_->align_max_shift_;
      double cur_pair_n_term_shift = cur_pair->getDiagonalHeader()
          ->getProtNTermShift();
      if (cur_pair_n_term_shift >= cur_pair_min_n_term_shift 
          && cur_pair_n_term_shift <= cur_pair_max_n_term_shift) {
        if (first_pair_->getSrc(s - 1) > shift_score) {
          shift_prev = first_pair_;
          shift_score = first_pair_->getSrc(s - 1);
        }
      }
    }

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
          || abs_shift < mng_->align_min_shift_
          || abs_shift > mng_->align_max_shift_) {
        continue;
      }
      double prev_score = prev_pair->getSrc(s - 1);
      if (abs_shift > mng_-> align_large_shift_thresh_) {
        prev_score = prev_score - mng_-> align_large_shift_panelty_;
      }
      if (prev_score > shift_score) {
        shift_prev = prev_pair;
        shift_score = prev_score;
      }
    }
  }
  return shift_prev;
}

void PSAlign::dp(SemiAlignTypePtr align_type) {
  dpPrep();
  for (unsigned int p = 1; p < dp_pairs_.size(); p++) {
    for (int s = 0; s <= mng_->n_unknown_shift_; s++) {
      DPPairPtr trunc_prev = getTruncPre(dp_pairs_[p], s, align_type);
      double trunc_score;
      if (trunc_prev == nullptr) {
        trunc_score = - std::numeric_limits<double>::max();
      } else {
        trunc_score = trunc_prev->getSrc(s);
      }

      DPPairPtr diag_prev = dp_pairs_[p]->getDiagPrev();
      double diag_score;
      if (diag_prev != nullptr) {
        diag_score = diag_prev->getSrc(s);
      } else {
        diag_score = - std::numeric_limits<double>::max();
      }
      DPPairPtr shift_prev = getShiftPre(dp_pairs_[p], p, s, align_type);
      double shift_score;
      if (shift_prev == nullptr) {
        shift_score =  - std::numeric_limits<double>::max();
      } else {
        shift_score = shift_prev->getSrc(s - 1);
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
    align_scores_.push_back(0.0);
    DiagonalHeaderPtrVec temp;
    backtrack_diagonals_.push_back(temp);
    backtrack_diagonals_[s] = backtrace(s);
  }
}

DiagonalHeaderPtrVec PSAlign::backtrace(int s) {
  DiagonalHeaderPtrVec list;
  DiagonalHeaderPtr cur_header;
  int cur_end = -1;
  int cur_bgn = -1;
  DPPairPtr p = last_pair_;

  align_scores_[s] = p->getSrc(s);
  if (p->getPre(s) == nullptr || p->getPre(s) == first_pair_) {
    return list;
  }

  while (p != first_pair_) {
    DPPairPtr pre = p->getPre(s);
    if (p == last_pair_) {
      cur_header = pre->getDiagonalHeader();
      cur_end = pre->getY() - 1;
    } else if (pre == first_pair_) {
      cur_bgn = p->getY();
      list.push_back(getDiagonalHeaderPtr(cur_header, cur_bgn, cur_end));
    } else {
      if (p->getType(s) == PATH_TYPE_SHIFT) {
        cur_bgn = p->getY();
        list.push_back(prot::getDiagonalHeaderPtr(cur_header, cur_bgn, cur_end));
        cur_header = pre->getDiagonalHeader();
        cur_end = pre->getY() - 1;
      }
    }
    if (p->getType(s) == PATH_TYPE_SHIFT) {
      s--;
    }
    p = pre;
  }
  DiagonalHeaderPtrVec sub_seg;
  for (unsigned int i = 0; i < list.size(); i++) {
    sub_seg.push_back(list[list.size() - 1 - i]);
  }
  return sub_seg;
}
} /* namespace prot */
