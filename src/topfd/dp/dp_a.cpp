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
#include "topfd/dp/vertex_util.hpp"
#include "topfd/dp/dp_a.hpp"

namespace toppic {

DpA::DpA(DeconvDataPtr data_ptr, MatchEnvPtr2D &win_envs_, 
         DpParaPtr dp_para_ptr, double score_error_tolerance):
    Dp(data_ptr, win_envs_, dp_para_ptr, score_error_tolerance) {
      initGraph();
      dp();
      backtrace();
    }

void DpA::initGraph() {
  // use win_num__ + 2 columns
  vertices_.resize(win_num_ + 2);
  MatchEnvPtrVec envs;
  vertices_[0] = vertex_util::getVertexAList(data_ptr_, -1, envs, envs, dp_para_ptr_);
  vertices_[1] = vertex_util::getVertexAList(data_ptr_, 0, envs, win_envs_[0], dp_para_ptr_);
  for (int i = 1; i < win_num_; i++) {
    vertices_[i + 1] = vertex_util::getVertexAList(data_ptr_, i, win_envs_[i - 1],
                                                  win_envs_[i], dp_para_ptr_);
  }
  vertices_[win_num_ + 1] = vertex_util::getVertexAList(data_ptr_, win_num_,
                                                       win_envs_[win_num_ - 1], envs, dp_para_ptr_);
  int cnt = 0;
  for (int i = 0; i <= win_num_ + 1; i++) {
    cnt += vertices_[i].size();
  }
  LOG_DEBUG("Vertex count " << cnt);
}

void DpA::dp() {
  int cnt = 0;
  for (int i = 1; i < win_num_ + 2; i++) {
    for (size_t j = 0; j < vertices_[i].size(); j++) {
      VertexAPtr cur_ver = vertices_[i][j];
      for (size_t k = 0; k < vertices_[i - 1].size(); k++) {
        VertexAPtr prev_ver = vertices_[i - 1][k];
        if (Vertex::checkConsist(prev_ver, cur_ver, dp_para_ptr_->max_env_num_per_peak_)) {
          cnt++;
          double new_score
              = Vertex::getShareScr(prev_ver, cur_ver, score_error_tolerance_);
          double cur_score
              = prev_ver->getScrA() + new_score;
          // LOG_DEBUG("i " << i << " j " << j << " k " << k << " new score " << new_score << " cur_score " << cur_score);
          if (cur_score > cur_ver->getScrA()) {
            cur_ver->setScrA(cur_score);
            cur_ver->setPreA(k);
          }
        }
      }
    }
  }
  LOG_DEBUG("Edge count :" << cnt);
}

// backtracking
void DpA::backtrace() {
  LOG_DEBUG("start backtrace ");
  int best_ver = -1;
  double best_score = - std::numeric_limits<double>::max();
  for (size_t i = 0; i < vertices_[win_num_ + 1].size(); i++) {
    double cur_score = vertices_[win_num_ + 1][i]->getScrA();
    if (cur_score > best_score) {
      best_ver = i;
      best_score = cur_score;
    }
  }
  // LOG_DEBUG("backtrace 1");
  for (int i = win_num_ + 1; i >= 1; i--) {
    // LOG_DEBUG("i " << i << " best ver " << best_ver);
    // LOG_DEBUG(" null " << (vertices_[i][best_ver]==nullptr));
    if (vertices_[i][best_ver]->getPreA() >= 0) {
      MatchEnvPtrVec prev_envs = vertices_[i][best_ver]->getPreMatchEnvs();
      addEnv(results_, prev_envs);
    }
    best_ver = vertices_[i][best_ver]->getPreA();
    if (best_ver < 0) {
      break;
    }
  }
  // LOG_DEBUG("backtrace 2");
  std::sort(results_.begin(), results_.end(), MatchEnv::cmpScoreDec);
}

}  // namespace toppic
