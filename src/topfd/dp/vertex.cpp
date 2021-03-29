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


#include "topfd/dp/vertex.hpp"

namespace toppic {

Vertex::Vertex(DpParaPtr dp_para_ptr, int bgn_peak, 
               int pre_win_peak_num, int cur_win_peak_num) {
  dp_para_ptr_ = dp_para_ptr;
  bgn_peak_ = bgn_peak;
  peak_num_ = pre_win_peak_num + cur_win_peak_num;
  cur_bgn_pos_ = pre_win_peak_num;
  peak_use_cnts_.resize(peak_num_, 0);

  pre_env_peak_pairs_.resize(peak_num_); 
  cur_env_peak_pairs_.resize(peak_num_); 
}

Vertex::Vertex(VertexPtr ptr) {
  dp_para_ptr_= ptr->dp_para_ptr_;
  bgn_peak_ = ptr->bgn_peak_;
  peak_num_ = ptr->peak_num_;
  cur_bgn_pos_ = ptr->cur_bgn_pos_;

  peak_use_cnts_ = ptr->peak_use_cnts_;

  prec_match_envs_ = ptr->prec_match_envs_;
  cur_match_envs_ = ptr->cur_match_envs_;

  pre_env_peak_pairs_ = ptr->pre_env_peak_pairs_;
  cur_env_peak_pairs_ = ptr->cur_env_peak_pairs_;
}

// check conexistance table 
bool Vertex::passDblIncrCheck(MatchEnvPtr env) {
  int id = env->getId();
  int exist_id;
  for (size_t i = 0; i < prec_match_envs_.size(); i++) {
    exist_id = prec_match_envs_[i]->getId();
    if (!dp_para_ptr_->coexist_table_[id][id - exist_id - 1]) {
      return false;
    }
  }

  for (size_t i = 0; i < cur_match_envs_.size(); i++) {
    exist_id = cur_match_envs_[i]->getId();
    if (!dp_para_ptr_->coexist_table_[id][id - exist_id - 1]) {
      return false;
    }
  }
  return true;
}

// check conexistance table 
bool Vertex::passDblIncrCheck(MatchEnvPtrVec &env_list) {
  for (size_t i = 0; i < env_list.size(); i++) {
    if (!passDblIncrCheck(env_list[i])) {
      return false;
    }
  }
  return true;
}

// Add an envelope to the previous window envelope list 
bool Vertex::addPreEnv(MatchEnvPtr env, int max_overlap) {
  if (dp_para_ptr_->check_double_increase_ && !passDblIncrCheck(env)) {
    return false;
  }
  RealEnvPtr real_env_ptr = env->getRealEnvPtr();
  //std::vector<int> idx_list = env->getRealEnvPtr()->getPeakIdxList();
  for (int i = 0; i < real_env_ptr->getPeakNum(); i++) {
    int peak_idx = real_env_ptr->getPeakIdx(i);
    if (real_env_ptr->isExist(i) && peak_idx >= bgn_peak_
        && peak_idx - bgn_peak_ < peak_num_) {

      peak_use_cnts_[peak_idx - bgn_peak_]++;
      if (peak_use_cnts_[peak_idx - bgn_peak_] > max_overlap) {
        return false;
      }
      EnvPeakPairPtr pair_ptr = std::make_shared<EnvPeakPair>(env, i);
      pre_env_peak_pairs_[peak_idx - bgn_peak_].push_back(pair_ptr);
    }
  }
  prec_match_envs_.push_back(env);
  return true;
}

// Add an envelope to the current window envelope list 
bool Vertex::addCurEnv(MatchEnvPtr env, int max_overlap) {
  if (dp_para_ptr_->check_double_increase_ && !passDblIncrCheck(env)) {
    return false;
  }
  RealEnvPtr real_env_ptr = env->getRealEnvPtr();
  for (int i = 0; i < real_env_ptr->getPeakNum(); i++) {
    int peak_idx = real_env_ptr->getPeakIdx(i);
    if (real_env_ptr->isExist(i) && peak_idx >= bgn_peak_
        && peak_idx - bgn_peak_ < peak_num_) {
      peak_use_cnts_[peak_idx - bgn_peak_]++;
      if (peak_use_cnts_[peak_idx - bgn_peak_] > max_overlap) {
        return false;
      }
      EnvPeakPairPtr pair_ptr = std::make_shared<EnvPeakPair>(env, i);
      cur_env_peak_pairs_[peak_idx - bgn_peak_].push_back(pair_ptr);
    }
  }
  cur_match_envs_.push_back(env);
  return true;
}

void Vertex::trim() {
  prec_match_envs_.shrink_to_fit();
  cur_match_envs_.shrink_to_fit();

  for (int i = 0; i < peak_num_; i++) {
    pre_env_peak_pairs_[i].shrink_to_fit();
    cur_env_peak_pairs_[i].shrink_to_fit();
  }
}

// check if two vertices in neighboring windows are consistent 
bool Vertex::checkConsist(VertexPtr pre, VertexPtr cur, int max_env_per_peak) {
  // if same size 
  if (pre->cur_match_envs_.size() != cur->prec_match_envs_.size()) {
    return false;
  }
  // if same envelope list 
  for (size_t i = 0; i < pre->cur_match_envs_.size(); i++) {
    if (pre->cur_match_envs_[i] != cur->prec_match_envs_[i]) {
      return false;
    }
  }
  // if no more than max overlap 
  for (size_t i = 0; i < pre->pre_env_peak_pairs_.size(); i++) {
    if (pre->pre_env_peak_pairs_[i].size() > 0) {
      int peak_idx = i + pre->bgn_peak_ - cur->bgn_peak_;
      if (peak_idx >= 0
          && (int)cur->peak_use_cnts_[peak_idx] + (int)pre->pre_env_peak_pairs_[i].size() 
          > max_env_per_peak) {
        return false;
      }
    }
  }
  if (pre->dp_para_ptr_->check_double_increase_) {
    MatchEnvPtrVec envs = cur->getCurMatchEnvs();
    if (!pre->passDblIncrCheck(envs)) {
      return false;
    }
  }
  return true;
}

// compute the score for peak sharing model 
double Vertex::getShareScr(VertexPtr pre, VertexPtr cur, double error_tole) {
  double score = 0;
  for (int cur_pnt = 0; cur_pnt < cur->cur_bgn_pos_; cur_pnt++) {
    // merge three peak pair list 
    int pre_pnt = cur_pnt + cur->bgn_peak_ - pre->bgn_peak_;
    int env_num = pre->pre_env_peak_pairs_[pre_pnt].size()
        + cur->pre_env_peak_pairs_[cur_pnt].size()
        + cur->cur_env_peak_pairs_[cur_pnt].size();
    EnvPeakPairPtrVec pairs(env_num);
    int cnt = 0;
    for (size_t i = 0; i < pre->pre_env_peak_pairs_[pre_pnt].size(); i++) {
      pairs[cnt] = pre->pre_env_peak_pairs_[pre_pnt][i];
      cnt++;
    }
    for (size_t i = 0; i < cur->pre_env_peak_pairs_[cur_pnt].size(); i++) {
      pairs[cnt] = cur->pre_env_peak_pairs_[cur_pnt][i];
      cnt++;
    }
    for (size_t i = 0; i < cur->cur_env_peak_pairs_[cur_pnt].size(); i++) {
      pairs[cnt] = cur->cur_env_peak_pairs_[cur_pnt][i];
      cnt++;
    }
    // compute sum of intensity 
    double inte_sum = 0;
    for (int i = 0; i < env_num; i++) {
      inte_sum += pairs[i]->getTheoIntensity();
    }
    // compute score 
    for (int i = 0; i < env_num; i++) {
      score += pairs[i]->getPeakScore(inte_sum, error_tole);
    }
  }
  return score;
}

}
