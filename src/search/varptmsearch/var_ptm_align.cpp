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

#include "search/varptmsearch/var_ptm_align.hpp"

#include <algorithm>
#include <limits>

#include "common/base/ptm.hpp"
#include "common/util/logger.hpp"
#include "ms/factory/extend_ms_factory.hpp"
#include "search/diag/diag_header_util.hpp"
#include "search/diag/diagonal_util.hpp"
#include "seq/alter.hpp"
#include "seq/mass_shift.hpp"
#include "seq/proteoform_factory.hpp"

namespace toppic {

namespace {

// whether one of the modifications modifies the residue
bool isModifiable(const ResiduePtr& res_ptr, const ModPtrVec& mod_ptrs) {
  for (size_t i = 0; i < mod_ptrs.size(); i++) {
    if (res_ptr->isSame(mod_ptrs[i]->getOriResiduePtr())) {
      return true;
    }
  }
  return false;
}

}  // namespace

VarPtmAlign::VarPtmAlign(const std::vector<double>& ms_masses,
                         const std::vector<double>& seq_masses,
                         const ResSeqPtr& sub_res_seq_ptr,
                         const DiagonalPtrVec& diagonal_ptrs,
                         const VarPtmSearchMngPtr& mng_ptr, int max_ptm_num)
    : ms_masses_(ms_masses),
      seq_masses_(seq_masses),
      sub_res_seq_ptr_(sub_res_seq_ptr),
      diagonal_ptrs_(diagonal_ptrs),
      mng_ptr_(mng_ptr),
      max_ptm_num_(max_ptm_num) {
  initMatchTable();
  initAllowShiftTable();
  initScoreTable();
}

void VarPtmAlign::initMatchTable() {
  for (size_t i = 0; i < diagonal_ptrs_.size(); i++) {
    std::vector<int> matches(seq_masses_.size(), 0);
    DiagPairPtrVec diag_pair_ptrs = diagonal_ptrs_[i]->getDiagPairPtrVec();
    // LOG_ERROR(i << " match number " << diag_pair_ptrs.size());
    for (size_t j = 0; j < diag_pair_ptrs.size(); j++) {
      int pos = diag_pair_ptrs[j]->getY();
      // LOG_ERROR(i << " " << j << " pos " << pos);
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
      ResiduePtr seq_res_ptr = sub_res_seq_ptr_->getResiduePtr(i - 1);
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
        scores_3d_[d][pos][ptm] =
            scores_3d_[d][pos - 1][ptm] + matches_2d_[d][pos];
      }
      for (ptm = 1; ptm <= var_ptm_num; ptm++) {
        int best_prev_d = d;
        int best_score = scores_3d_[d][pos - 1][ptm];
        for (size_t k = 0; k < mng_ptr_->diag_prev_idxes_[d].size(); k++) {
          int shift_idx = mng_ptr_->diag_prev_shift_idxes_[d][k];
          bool allow_shift = allow_shift_2d_[pos][shift_idx];
          if (allow_shift) {
            int prev_d = mng_ptr_->diag_prev_idxes_[d][k];
            int score = scores_3d_[prev_d][pos - 1][ptm - 1];
            if (score > best_score) {
              best_score = score;
              best_prev_d = prev_d;
            }
          }
        }
        // add match score for the current position
        best_score += matches_2d_[d][pos];
        prevs_3d_[d][pos][ptm] = best_prev_d;
        scores_3d_[d][pos][ptm] = best_score;
        /*
        if (pos == seq_masses_.size() - 1 && best_score > 0) {
          LOG_ERROR("d " << d << " ptm " << ptm << " best score " <<
        best_score);
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
  int var_ptm_num = std::min(max_ptm_num_, mng_ptr_->var_ptm_num_);
  int pos = seq_masses_.size() - 1;
  for (size_t i = 0; i < diagonal_ptrs_.size(); i++) {
    LOG_DEBUG("diag " << i << " shift "
                      << diagonal_ptrs_[i]->getHeader()->getProtNTermShift()
                      << " c term match "
                      << diagonal_ptrs_[i]->getHeader()->isProtCTermMatch());
    // check if the diagonal is matched to the C-terminus
    if (diagonal_ptrs_[i]->getHeader()->isProtCTermMatch()) {
      for (int j = 0; j <= var_ptm_num; j++) {
        int score = scores_3d_[i][pos][j];
        if (score > best_score_) {
          best_score_ = score;
          best_d_ = i;
          best_ptm_num_ = j;
        }
      }
    }
  }
  LOG_DEBUG("best d " << best_d_ << " best ptm num " << best_ptm_num_
                      << " best score " << best_score_);
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
      list.push_back(
          diag_header_util::geneDiagHeaderPtr(cur_bgn, cur_end, cur_header));

      // find mod
      ModPtr cur_mod_ptr;
      int res_pos = cur_bgn - 1;
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
  list.push_back(
      diag_header_util::geneDiagHeaderPtr(cur_bgn, cur_end, cur_header));
  std::reverse(list.begin(), list.end());
  backtrack_diag_header_ptrs_ = list;
  std::reverse(backtrack_mod_ptrs_.begin(), backtrack_mod_ptrs_.end());
}

MassShiftPtrVec VarPtmAlign::geneShiftVec(
    const ProteoformPtr& sub_proteo_ptr, DiagHeaderPtrVec& ori_header_ptrs,
    DiagHeaderPtrVec& refined_header_ptrs) {
  DiagHeaderPtrVec filtered_header_ptrs;
  // get non-null headers
  for (size_t i = 0; i < refined_header_ptrs.size(); i++) {
    LOG_DEBUG("refine header " << i << " null "
                               << (refined_header_ptrs[i] == nullptr));
    if (refined_header_ptrs[i] != nullptr) {
      filtered_header_ptrs.push_back(refined_header_ptrs[i]);
    }
  }
  if (filtered_header_ptrs.size() == 0) {
    LOG_ERROR("Empty diagonal header");
    exit(EXIT_FAILURE);
  }
  int first_pos = filtered_header_ptrs[0]->getTruncFirstResPos();
  int last_pos = filtered_header_ptrs[filtered_header_ptrs.size() - 1]
                     ->getTruncLastResPos();
  LOG_DEBUG("first pos " << first_pos << " last pos " << last_pos);
  MassShiftPtrVec shift_ptrs = diag_header_util::getDiagonalMassChanges(
      filtered_header_ptrs, first_pos, last_pos, AlterType::VARIABLE);
  // add AlterPtr to shifts
  int idx = 0;
  AlterPtrVec alter_ptr_vec;
  for (size_t i = 1; i < ori_header_ptrs.size(); i++) {
    DiagHeaderPtr cur_ptr = ori_header_ptrs[i];
    DiagHeaderPtr pre_ptr = ori_header_ptrs[i - 1];
    double shift = cur_ptr->getProtNTermShift() - pre_ptr->getProtNTermShift();
    int left_bp_pos = shift_ptrs[idx]->getLeftBpPos();
    int right_bp_pos = shift_ptrs[idx]->getRightBpPos();
    ModPtr mod_ptr = backtrack_mod_ptrs_[i - 1];
    // find modifiable residues: those of every modification of the PTM
    // (Deamidated on N and on Q), not only of the residue the alignment chose
    const ModPtrVec& ptm_mod_ptrs = mng_ptr_->getPtmModPtrVec(mod_ptr);
    while (left_bp_pos < right_bp_pos) {
      if (isModifiable(sub_res_seq_ptr_->getResiduePtr(left_bp_pos),
                       ptm_mod_ptrs)) {
        break;
      }
      left_bp_pos++;
    }
    while (right_bp_pos > left_bp_pos) {
      if (isModifiable(sub_res_seq_ptr_->getResiduePtr(right_bp_pos - 1),
                       ptm_mod_ptrs)) {
        break;
      }
      right_bp_pos--;
    }

    AlterPtr a = std::make_shared<Alter>(left_bp_pos, right_bp_pos,
                                         AlterType::VARIABLE, shift, mod_ptr);
    alter_ptr_vec.push_back(a);

    if (refined_header_ptrs[i] != nullptr) {
      shift_ptrs[idx]->setAlterPtrVec(alter_ptr_vec);
      // the shift covers the residues that can carry its PTMs, not the whole
      // gap between the refined diagonals: the cleavages outside them keep
      // their ions
      shift_ptrs[idx]->narrowRangeToAlters();
      alter_ptr_vec.clear();
      idx++;
    }
  }
  // a shift with as many candidate residues as PTMs becomes one shift per
  // PTM on its residue
  MassShiftPtrVec result_ptrs;
  for (size_t i = 0; i < shift_ptrs.size(); i++) {
    MassShiftPtrVec cur_ptrs = splitDeterminedShift(shift_ptrs[i]);
    result_ptrs.insert(result_ptrs.end(), cur_ptrs.begin(), cur_ptrs.end());
  }
  return result_ptrs;
}

MassShiftPtrVec VarPtmAlign::splitDeterminedShift(
    const MassShiftPtr& shift_ptr) {
  MassShiftPtrVec result_ptrs{shift_ptr};
  const AlterPtrVec& alter_ptrs = shift_ptr->getAlterPtrVec();
  if (alter_ptrs.size() < 2) {
    return result_ptrs;
  }
  // the alterations of each PTM and the candidate residues of the PTM
  std::vector<AlterPtrVec> ptm_alter_ptrs;
  std::vector<std::vector<int>> ptm_positions;
  for (size_t i = 0; i < alter_ptrs.size(); i++) {
    const PtmPtr& ptm_ptr =
        alter_ptrs[i]->getModPtr()->getModResiduePtr()->getPtmPtr();
    size_t g = 0;
    while (g < ptm_alter_ptrs.size() &&
           !ptm_alter_ptrs[g][0]->getModPtr()->getModResiduePtr()
                ->getPtmPtr()->isSame(ptm_ptr)) {
      g++;
    }
    if (g == ptm_alter_ptrs.size()) {
      const ModPtrVec& ptm_mod_ptrs =
          mng_ptr_->getPtmModPtrVec(alter_ptrs[i]->getModPtr());
      std::vector<int> positions;
      for (int pos = shift_ptr->getLeftBpPos();
           pos < shift_ptr->getRightBpPos(); pos++) {
        if (isModifiable(sub_res_seq_ptr_->getResiduePtr(pos), ptm_mod_ptrs)) {
          positions.push_back(pos);
        }
      }
      ptm_alter_ptrs.push_back(AlterPtrVec());
      ptm_positions.push_back(positions);
    }
    ptm_alter_ptrs[g].push_back(alter_ptrs[i]);
  }
  // determined: as many candidate residues as alterations for each PTM, and
  // no residue a candidate of two PTMs
  std::vector<int> used_positions;
  for (size_t g = 0; g < ptm_alter_ptrs.size(); g++) {
    if (ptm_positions[g].size() != ptm_alter_ptrs[g].size()) {
      return result_ptrs;
    }
    for (size_t j = 0; j < ptm_positions[g].size(); j++) {
      if (std::find(used_positions.begin(), used_positions.end(),
                    ptm_positions[g][j]) != used_positions.end()) {
        return result_ptrs;
      }
      used_positions.push_back(ptm_positions[g][j]);
    }
  }
  // one shift per alteration, with the modification of its residue
  // (Deamidated on N or on Q)
  result_ptrs.clear();
  for (size_t g = 0; g < ptm_alter_ptrs.size(); g++) {
    const ModPtrVec& ptm_mod_ptrs =
        mng_ptr_->getPtmModPtrVec(ptm_alter_ptrs[g][0]->getModPtr());
    for (size_t j = 0; j < ptm_alter_ptrs[g].size(); j++) {
      int pos = ptm_positions[g][j];
      const ResiduePtr& res_ptr = sub_res_seq_ptr_->getResiduePtr(pos);
      ModPtr mod_ptr = ptm_alter_ptrs[g][j]->getModPtr();
      for (size_t k = 0; k < ptm_mod_ptrs.size(); k++) {
        if (res_ptr->isSame(ptm_mod_ptrs[k]->getOriResiduePtr())) {
          mod_ptr = ptm_mod_ptrs[k];
          break;
        }
      }
      AlterPtr alter_ptr = std::make_shared<Alter>(
          pos, pos + 1, AlterType::VARIABLE, ptm_alter_ptrs[g][j]->getMass(),
          mod_ptr);
      result_ptrs.push_back(std::make_shared<MassShift>(alter_ptr));
    }
  }
  std::sort(result_ptrs.begin(), result_ptrs.end(), MassShift::cmpPosInc);
  return result_ptrs;
}

PrsmPtr VarPtmAlign::geneResult(const ProteoformPtr& sub_proteo_ptr,
                                DeconvMsPtrVec& deconv_ms_ptr_vec,
                                const PrsmParaPtr& prsm_para_ptr) {
  if (backtrack_diag_header_ptrs_.size() == 0) {
    return nullptr;
  }

  // the precursor mass of the proteoform with the variable PTMs of the
  // C-terminal diagonal; the backtracked headers take their C-terminal
  // shifts from it before their ranges are refined
  double refine_prec_mass = diagonal_util::refineVarPtmPrecursorAndHeaderShift(
      sub_proteo_ptr, backtrack_diag_header_ptrs_,
      mng_ptr_->shift_list_[best_d_]);

  SpParaPtr sp_para_ptr = prsm_para_ptr->getSpParaPtr();
  ExtendMsPtrVec refine_ms_ptr_vec = extend_ms_factory::geneMsThreePtrVec(
      deconv_ms_ptr_vec, sp_para_ptr, refine_prec_mass);

  double min_mass = prsm_para_ptr->getSpParaPtr()->getMinMass();
  DiagHeaderPtrVec refined_header_ptrs =
      diagonal_util::refineVarPtmHeadersBgnEnd(
          sub_proteo_ptr, refine_ms_ptr_vec, backtrack_diag_header_ptrs_,
          min_mass);
  if (refined_header_ptrs.size() == 0) {
    return nullptr;
  }

  MassShiftPtrVec shift_vec = geneShiftVec(
      sub_proteo_ptr, backtrack_diag_header_ptrs_, refined_header_ptrs);

  sub_proteo_ptr->addMassShiftPtrVec(shift_vec);

  return std::make_shared<Prsm>(sub_proteo_ptr, deconv_ms_ptr_vec,
                                refine_prec_mass,
                                prsm_para_ptr->getSpParaPtr());
}

} /* namespace toppic */
