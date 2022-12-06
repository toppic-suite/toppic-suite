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
                         ResSeqPtr sub_res_seq_ptr,
                         const DiagonalPtrVec &diagonal_ptrs,
                         VarPtmSearchMngPtr mng_ptr):
  ms_masses_(ms_masses),
  seq_masses_(seq_masses),
  sub_res_seq_ptr_(sub_res_seq_ptr),
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
    //LOG_ERROR(i << " match number " << diag_pair_ptrs.size());
    for (size_t j = 0; j < diag_pair_ptrs.size(); j++) {
      int pos = diag_pair_ptrs[j]->getY();
      //LOG_ERROR(i << " " << j << " pos " << pos);
      matches[pos] = 1;
    }
    matches_2d_.push_back(matches);
  }
}

void VarPtmAlign::initAllowShiftTable() {
  size_t shift_num = mng_ptr_->single_shift_list_.size();
  for (size_t i = 0; i < seq_masses_.size(); i++) {
    std::vector<bool> allow_shifts(shift_num, false);
    if (i >= 1) {
      ResiduePtr seq_res_ptr = sub_res_seq_ptr_->getResiduePtr(i-1);
      for (size_t j = 0; j < shift_num; j++) {
        //LOG_DEBUG("shift " << j << " mod num " << mng_ptr_->res_ptr_vec_2d_[j].size());
        for (size_t k = 0; k < mng_ptr_->res_ptr_vec_2d_[j].size(); k++) {
          ResiduePtr mod_res_ptr = mng_ptr_->res_ptr_vec_2d_[j][k];
          if (seq_res_ptr->isSame(mod_res_ptr)) {
            allow_shifts[j] = true;
            break;
          }
        }
      }
    }
    //LOG_ERROR("Allow shift " << i << " " << allow_shifts[0] << " " << allow_shifts[1]);
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
  backtrack();
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
        /*
        if (pos == seq_masses_.size() - 1 && best_score > 0) {
          LOG_ERROR("d " << d << " ptm " << ptm << " best score " << best_score);
        }
        */
      }
    }
  }
}

void VarPtmAlign::backtrack() {
  best_score_ = -1;
  best_d_ = 0;
  best_ptm_num_ = 0;
  size_t var_ptm_num = mng_ptr_->var_ptm_num_;
  int pos = seq_masses_.size() -1 ;
  for (size_t i = 0; i < diagonal_ptrs_.size(); i++) {
    LOG_DEBUG("diag " << i << " shift " << diagonal_ptrs_[i]->getHeader()->getProtNTermShift()
              << " c term match " << diagonal_ptrs_[i]->getHeader()->isProtCTermMatch())
    // check if the diagonal is matched to the C-terminus
    if (diagonal_ptrs_[i]->getHeader()->isProtCTermMatch()) {
      for (size_t j = 0; j <= var_ptm_num; j++) {
        int score = scores_3d_[i][pos][j];
        if (score > best_score_) {
          best_score_ = score;
          best_d_ = i;
          best_ptm_num_ = j;
        }
      }
    }
  }
  LOG_DEBUG("best d " << best_d_ << " best ptm num " << best_ptm_num_ << " best score " << best_score_);
  if (best_score_ > 0) {
    backtrack(best_d_, best_ptm_num_);
  }
}

void VarPtmAlign::backtrack(int d, int ptm) {

  DiagHeaderPtrVec list;
  int cur_d = d;
  int cur_pos = seq_masses_.size() - 1;
  int cur_ptm = ptm;
  if (scores_3d_[cur_d][cur_pos][cur_ptm] <= 0) {
    return; 
  }

  int cur_end = cur_pos;
  int cur_bgn = cur_pos;
  DiagHeaderPtr cur_header;

  // LOG_DEBUG("backtrace start ");
  while (cur_pos != 0) {
    int pre_d = prevs_3d_[cur_d][cur_pos][cur_ptm];
    int pre_pos = cur_pos - 1;
    int pre_ptm = cur_ptm;
    if (pre_d != cur_d) {
      pre_ptm = cur_ptm - 1;
    }

    if (pre_d != cur_d) {
      cur_bgn = cur_pos;
      cur_header = diagonal_ptrs_[cur_d]->getHeader();
      list.push_back(diag_header_util::geneDiagHeaderPtr(cur_bgn, cur_end, cur_header));

      // find mod
      ModPtr cur_mod_ptr;
      int res_pos = cur_bgn  - 1;
      ResiduePtr seq_res_ptr = sub_res_seq_ptr_->getResiduePtr(res_pos);
      int shift_idx = mng_ptr_->diag_matrix_shift_idxes_[cur_d][pre_d];

      ModPtrVec cand_mod_ptrs = mng_ptr_->mod_ptr_vec_2d_[shift_idx];
      for (size_t k = 0; k < cand_mod_ptrs.size(); k++) {
        ResiduePtr mod_res_ptr = cand_mod_ptrs[k]->getOriResiduePtr();
        if (seq_res_ptr->isSame(mod_res_ptr)) {
          cur_mod_ptr = cand_mod_ptrs[k];
          break;
        }
      }
      if (cur_mod_ptr == nullptr) {
        LOG_ERROR("No matched modification is found!");
        exit(EXIT_FAILURE);
      }
      backtrack_mod_ptrs_.push_back(cur_mod_ptr);

      cur_end = pre_pos;
      cur_bgn = pre_pos;
    }
    cur_d = pre_d;
    cur_pos = pre_pos;
    cur_ptm = pre_ptm;
  }
  // first diagonal;
  cur_bgn = cur_pos;
  cur_header = diagonal_ptrs_[cur_d]->getHeader();
  list.push_back(diag_header_util::geneDiagHeaderPtr(cur_bgn, cur_end, cur_header));
  std::reverse(list.begin(), list.end());
  backtrack_diag_header_ptrs_ = list;
  std::reverse(backtrack_mod_ptrs_.begin(), backtrack_mod_ptrs_.end());
}

MassShiftPtrVec VarPtmAlign::geneShiftVec(ProteoformPtr sub_proteo_ptr,
                                          DiagHeaderPtrVec &ori_header_ptrs,
                                          DiagHeaderPtrVec &refined_header_ptrs) {
  DiagHeaderPtrVec filtered_header_ptrs;
  // get non-null headers
  for (size_t i = 0; i < refined_header_ptrs.size(); i++) {
    LOG_DEBUG("refine header " << i << " null " << (refined_header_ptrs[i] == nullptr));
    if (refined_header_ptrs[i] != nullptr) {
      filtered_header_ptrs.push_back(refined_header_ptrs[i]);
    }
  }
  if (filtered_header_ptrs.size() == 0) {
    LOG_ERROR("Empty diagonal header");
    exit(EXIT_FAILURE);
  }
  int first_pos = filtered_header_ptrs[0]->getTruncFirstResPos();
  int last_pos = filtered_header_ptrs[filtered_header_ptrs.size()-1]->getTruncLastResPos();
  LOG_DEBUG("first pos " << first_pos << " last pos " << last_pos); 
  MassShiftPtrVec shift_ptrs = diag_header_util::getDiagonalMassChanges(filtered_header_ptrs, first_pos,
                                                                        last_pos, AlterType::VARIABLE);
  // add AlterPtr to shifts
  int idx = 0;
  AlterPtrVec alter_ptr_vec;
  for (size_t i = 1; i < ori_header_ptrs.size(); i++) {
    DiagHeaderPtr cur_ptr = ori_header_ptrs[i];
    DiagHeaderPtr pre_ptr = ori_header_ptrs[i-1];
    double shift = cur_ptr->getProtNTermShift() - pre_ptr->getProtNTermShift(); 
    int left_bp_pos = shift_ptrs[idx]->getLeftBpPos();
    int right_bp_pos = shift_ptrs[idx]->getRightBpPos(); 
    ModPtr mod_ptr = backtrack_mod_ptrs_[i-1];
    ResiduePtr mod_res_ptr = mod_ptr->getOriResiduePtr();
    // find modifiable residues
    while (left_bp_pos < right_bp_pos) {
      ResiduePtr seq_res_ptr = sub_res_seq_ptr_->getResiduePtr(left_bp_pos);
      if (seq_res_ptr->isSame(mod_res_ptr)) {
        break;
      }
      left_bp_pos ++;
    }
    while (right_bp_pos > left_bp_pos) {
      ResiduePtr seq_res_ptr = sub_res_seq_ptr_->getResiduePtr(right_bp_pos - 1);
      if (seq_res_ptr->isSame(mod_res_ptr)) {
        break;
      }
      right_bp_pos --;
    }

    AlterPtr a = std::make_shared<Alter>(left_bp_pos, 
                                         right_bp_pos,
                                         AlterType::VARIABLE,
                                         shift, 
                                         mod_ptr);
    alter_ptr_vec.push_back(a);

    if (refined_header_ptrs[i] != nullptr) {
      shift_ptrs[idx]->setAlterPtrVec(alter_ptr_vec);
      alter_ptr_vec.clear();
      idx ++;
    }
  }
  return shift_ptrs;
}



PrsmPtr VarPtmAlign::geneResult(ProteoformPtr sub_proteo_ptr,
                                DeconvMsPtrVec &deconv_ms_ptr_vec,
                                PrsmParaPtr prsm_para_ptr) {
  if (backtrack_diag_header_ptrs_.size() == 0) {
    return nullptr;
  }

  double refine_prec_mass = sub_proteo_ptr->getMass() + mng_ptr_->shift_list_[best_d_];

  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  ExtendMsPtrVec refine_ms_ptr_vec = extend_ms_factory::geneMsThreePtrVec(deconv_ms_ptr_vec,
                                                                          sp_para_ptr, 
                                                                          refine_prec_mass);

  double min_mass = prsm_para_ptr->getSpParaPtr()->getMinMass();
  DiagHeaderPtrVec refined_header_ptrs 
    = diagonal_util::refineVarPtmHeadersBgnEnd(sub_proteo_ptr, 
                                               refine_ms_ptr_vec,
                                               backtrack_diag_header_ptrs_, 
                                               min_mass);
  if (refined_header_ptrs.size() == 0) {
    return nullptr;
  }

  MassShiftPtrVec shift_vec = geneShiftVec(sub_proteo_ptr,
                                           backtrack_diag_header_ptrs_, 
                                           refined_header_ptrs);

  sub_proteo_ptr->addMassShiftPtrVec(shift_vec);

  return std::make_shared<Prsm>(sub_proteo_ptr, deconv_ms_ptr_vec, refine_prec_mass,
                                prsm_para_ptr->getSpParaPtr());
}

} /* namespace toppic */
