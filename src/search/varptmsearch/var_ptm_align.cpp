//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#include <limits>
#include <algorithm>

#include "common/util/logger.hpp"
#include "seq/proteoform_factory.hpp"
#include "ms/factory/extend_ms_factory.hpp"
#include "search/diag/diagonal_util.hpp"
#include "search/diag/diag_header_util.hpp"
#include "search/varptmsearch/var_ptm_align.hpp"

namespace toppic {

VarPtmAlign::VarPtmAlign(const std::vector<double> &ms_masses,
                         const std::vector<double> &seq_masses,
                         ResSeqPtr res_seq_ptr,
                         const DiagonalPtrVec &diagonal_ptrs,
                         VarPtmSearchMngPtr mng_ptr):
  ms_masses_(ms_masses),
  seq_masses_(seq_masses),
  res_seq_ptr_(res_seq_ptr),
  diagonal_ptrs_(diagonal_ptrs),
  mng_ptr_(mng_ptr) {
    initMatchTable();
    initAllowShiftTable();
    initScoreTable();
  }

void VarPtmAlign::initMatchTable() {
  for (size_t i = 0; i < diagonal_ptrs_.size(); i++) {
    std::vector<int> matches(seq_masses_.size(), 0);
    DiagPairPtrVec diag_pair_ptrs = diagonal_ptrs_[i]->getDiagPairPtrVec();
    for (size_t j = 0; j < diag_pair_ptrs.size(); j++) {
      int pos = diag_pair_ptrs[j]->getY();
      LOG_DEBUG("pos" << pos);
      matches[j] = 1;
    }
    matches_2d_.push_back(matches);
  }
}

void VarPtmAlign::initAllowShiftTable() {
  size_t shift_num = mng_ptr_->single_shift_list_.size();
  for (size_t i = 0; i < seq_masses_.size(); i++) {
    std::vector<bool> allow_shifts(shift_num, false);
    if (i >= 1) {
      ResiduePtr seq_res_ptr = res_seq_ptr_->getResiduePtr(i-1);
      for (size_t j = 0; j < shift_num; j++) {
        for (size_t k = 0; k < mng_ptr_->res_ptr_vec_2d_[j].size(); k++) {
          ResiduePtr mod_res_ptr = mng_ptr_->res_ptr_vec_2d_[j][k];
          if (seq_res_ptr->isSame(mod_res_ptr)) {
            allow_shifts[j] = true;
            break;
          }
        }
      }
    }
    allow_shift_2d_.push_back(allow_shifts); 
  }
}

void VarPtmAlign::initScoreTable() {
  int min_int = std::numeric_limits<int>::min();
  int var_ptm_num = mng_ptr_->var_ptm_num_;
  for (size_t i = 0; i < diagonal_ptrs_.size(); i++) {
    std::vector<std::vector<int>> row_scores;  
    std::vector<std::vector<int>> row_prevs;  
    for (size_t j = 0; j < seq_masses_.size(); j++) {
      std::vector<int> cell_scores(var_ptm_num + 1, min_int); 
      std::vector<int> cell_prevs(var_ptm_num + 1, -1); 
      row_scores.push_back(cell_scores);
      row_prevs.push_back(cell_prevs);
    }
    scores_3d_.push_back(row_scores); 
    prevs_3d_.push_back(row_prevs);
  }
}


void VarPtmAlign::compute() {
  dp();
  //backtrace();
}

void VarPtmAlign::dp() {
  scores_3d_[0][0][0] = 0;
  int var_ptm_num = mng_ptr_->var_ptm_num_;
  for (size_t pos = 1; pos < seq_masses_.size(); pos++) {
    for (size_t d = 0; d < diagonal_ptrs_.size(); d++) {
      int ptm = 0;
      // compute for zero ptm
      if (d == 0) {
        prevs_3d_[d][pos][ptm] = d; 
        scores_3d_[d][pos][ptm] = scores_3d_[d][pos-1][ptm] + matches_2d_[d][pos];
      }
      for (ptm = 1; ptm <= var_ptm_num; ptm++) {
        int best_prev_d = d;
        int best_score = scores_3d_[d][pos-1][ptm];
        for (size_t k = 0; k < mng_ptr_->diag_prev_idxes_[d].size(); k++) {
          int shift_idx = mng_ptr_->diag_prev_shift_idxes_[d][k];
          bool allow_shift = allow_shift_2d_[pos][shift_idx];
          if (allow_shift) {
            int prev_d = mng_ptr_->diag_prev_idxes_[d][k];
            int score = scores_3d_[prev_d][pos-1][ptm-1];
            if (score > best_score) {
              best_score = score;
              best_prev_d = prev_d;
            }
          }
        }
        // add match score for the current position
        best_score +=  matches_2d_[d][pos];
        prevs_3d_[d][pos][ptm] = best_prev_d;
        scores_3d_[d][pos][ptm] = best_score;
      }
    }
  }
}

DiagHeaderPtrVec VarPtmAlign::backtrack(int d, int s) {
  DiagHeaderPtrVec list;
  int cur_d = d;
  int cur_pos = seq_masses_.size();
  int cur_ptm = s
  if (scores_3d_[cur_d][cur_pos][cur_s] <= 0) {
    return list;
  }

  DiagHeaderPtr cur_header;
  int cur_end = cur_pos;
  int cur_bgn = cur_pos;

  // LOG_DEBUG("backtrace start ");
  while (cur_pos != 0) {
    int pre_d = prevs_3d_[cur_d][cur_pos][cur_ptm];
    int pre_pos = cur_pos - 1;
    int pre_ptm = cur_ptm;
    if (pre_d != cur_d) {
      pre_ptm = cur_ptm - 1;
    }
    if (p == last_pair_ptr_) {
      cur_header = pre->getDiagHeader();
      cur_end = pre->getY();
    } else if (pre == first_pair_ptr_) {
      cur_bgn = p->getY();
      list.push_back(diag_header_util::geneDiagHeaderPtr(cur_bgn, cur_end, cur_header));
    } else {
      if (p->getType(s) == path_type::TYPE_SHIFT) {
        cur_bgn = p->getY();
        list.push_back(diag_header_util::geneDiagHeaderPtr(cur_bgn, cur_end, cur_header));
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

/*
inline DpPairPtr PsAlign::getTruncPre(DpPairPtr cur_pair_ptr, int s,
                                      ProteoformTypePtr align_type_ptr) {
  DpPairPtr trunc_prev_ptr;
  if (cur_pair_ptr == last_pair_ptr_) {
    double trunc_score = - std::numeric_limits<double>::max();
    for (size_t i = 0; i < segment_end_pair_ptrs_.size(); i++) {
      DpPairPtr prev_pair_ptr = segment_end_pair_ptrs_[i];
      if (align_type_ptr == ProteoformType::COMPLETE || align_type_ptr == ProteoformType::SUFFIX) {
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
      if (align_type_ptr == ProteoformType::COMPLETE
          || align_type_ptr == ProteoformType::PREFIX) {
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

DpPairPtr PsAlign::getShiftPre(int p, int s) {
  if (s < 1) {
    return nullptr;
  }
  // last pair can not shift, only trunc is allowed;
  DpPairPtr shift_prev = nullptr;
  double shift_score = -std::numeric_limits<double>::max();

  // cur pair is not last pair
  // first pair can not shift, so q starts from 1
  for (size_t d = 0; d < idxes_[p].size(); d++) {
    if (idxes_[p][d] >= 0) {
      DpPairPtr prev_pair_ptr = dp_pair_ptrs_[idxes_[p][d]];
      double prev_score = prev_pair_ptr->getScore(s - 1);
      if (penalties_[p][d]) {
        prev_score = prev_score - para_ptr_-> align_large_shift_panelty_;
      }
      if (prev_score > shift_score) {
        shift_prev = prev_pair_ptr;
        shift_score = prev_score;
      }
    }
  }
  return shift_prev;
}

PrsmPtr PsAlign::geneResult(int shift_num, ProteoformPtr proteo_ptr,
                            DeconvMsPtrVec &deconv_ms_ptr_vec,
                            ExtendMsPtrVec &ms_three_ptr_vec,
                            PrsmParaPtr prsm_para_ptr) {
  DiagHeaderPtrVec header_ptrs = getDiagHeaders(shift_num);
  if (header_ptrs.size() == 0) {
    return nullptr;
  }
  int first_pos = header_ptrs[0]->getTruncFirstResPos();
  int last_pos = header_ptrs[header_ptrs.size()-1]->getTruncLastResPos();
  ProteoformPtr sub_proteo_ptr = proteoform_factory::geneSubProteoform(proteo_ptr, 
                                                                       proteo_ptr->getFastaSeqPtr(),
                                                                       first_pos, last_pos);

  double min_mass = prsm_para_ptr->getSpParaPtr()->getMinMass();
  double ppo = prsm_para_ptr->getSpParaPtr()->getPeakTolerancePtr()->getPpo();

  double refine_prec_mass 
    = diagonal_util::refinePrecursorAndHeaderShift(proteo_ptr, ms_three_ptr_vec,
                                                   header_ptrs, ppo, min_mass,
                                                   para_ptr_->refine_prec_step_width_);

  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  ExtendMsPtrVec refine_ms_ptr_vec = extend_ms_factory::geneMsThreePtrVec(deconv_ms_ptr_vec,
                                                                          sp_para_ptr, refine_prec_mass);

  DiagHeaderPtrVec refined_header_ptrs = diagonal_util::refineHeadersBgnEnd(proteo_ptr, 
                                                                            refine_ms_ptr_vec,
                                                                            header_ptrs, 
                                                                            min_mass);

  if (refined_header_ptrs.size() == 0) {
    return nullptr;
  }

  MassShiftPtrVec shifts = diag_header_util::getDiagonalMassChanges(refined_header_ptrs, first_pos,
                                                                    last_pos, AlterType::UNEXPECTED);
  sub_proteo_ptr->addMassShiftPtrVec(shifts);

  return std::make_shared<Prsm>(sub_proteo_ptr, deconv_ms_ptr_vec, refine_prec_mass,
                                prsm_para_ptr->getSpParaPtr());
}
*/

} /* namespace toppic */
