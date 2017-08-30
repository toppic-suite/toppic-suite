// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include <limits>
#include <cstddef>
#include <algorithm>

#include "base/logger.hpp"
#include "feature/vertex_base.hpp"
#include "feature/dp_a.hpp"

namespace prot {

DpA::DpA(DeconvDataPtr data_ptr, MatchEnvPtr2D &win_envs_, FeatureMngPtr mng_ptr):
    Dp(data_ptr, win_envs_, mng_ptr) {
      initGraph();
      dp();
      backtrace();
    }

void DpA::initGraph() {
  // use win_num__ + 2 columns
  vertices_.resize(win_num_ + 2);
  MatchEnvPtrVec envs;
  vertices_[0] = VertexBase::getVertexAList(data_ptr_, -1, envs, envs, mng_ptr_);
  vertices_[1] = VertexBase::getVertexAList(data_ptr_, 0, envs, win_envs_[0], mng_ptr_);
  for (int i = 1; i < win_num_; i++) {
    vertices_[i + 1] = VertexBase::getVertexAList(data_ptr_, i, win_envs_[i - 1],
                                                  win_envs_[i], mng_ptr_);
  }
  vertices_[win_num_ + 1] = VertexBase::getVertexAList(data_ptr_, win_num_,
                                                       win_envs_[win_num_ - 1], envs, mng_ptr_);
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
        if (Vertex::checkConsist(prev_ver, cur_ver, mng_ptr_->max_env_num_per_peak_)) {
          cnt++;
          double new_score
              = Vertex::getShareScr(prev_ver, cur_ver, mng_ptr_->score_error_tolerance_);
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

}  // namespace prot
