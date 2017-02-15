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


#include "feature/vertex_base.hpp"

namespace prot {

int VertexBase::getBgnPeak(int pre_win, DeconvDataPtr data) {
  int bgn_peak;
  if (pre_win < 0) {
    bgn_peak = 0;
  } else {
    bgn_peak = data->getBgnPeak(pre_win);
  }
  return bgn_peak;
}

int VertexBase::getWinPkNum(int win, DeconvDataPtr data) {
  if (win < 0 || win >= data->getWinNum()) {
    return 0;
  } else {
    return data->getIntervalPeakNum(win);
  }
}

// add an empty vertexA 
void VertexBase::addEmptyVertexA(FeatureMngPtr mng_ptr, VertexAPtrVec &result, 
                                 DeconvDataPtr data, int cur_win) {
  int pre_win = cur_win - 1;
  int bgn_peak = getBgnPeak(pre_win, data);
  int pre_win_peak_num = getWinPkNum(pre_win, data);
  int cur_win_peak_num = getWinPkNum(cur_win, data);
  VertexAPtr ptr(new VertexA(mng_ptr, bgn_peak, pre_win_peak_num, cur_win_peak_num));
  ptr->trim();
  result.push_back(ptr);
}

VertexAPtrVec VertexBase::getVertexAList(DeconvDataPtr data, int cur_win, 
                                         MatchEnvPtrVec prev_match_envs, 
                                         MatchEnvPtrVec cur_match_envs, FeatureMngPtr mng_ptr) {
  VertexAPtrVec result;
  // add an empty vertex 
  addEmptyVertexA(mng_ptr, result, data, cur_win);
  int cur_size = 1;
  // add envelopes in previous window 
  for (size_t i = 0; i < prev_match_envs.size(); i++) {
    cur_size = result.size();
    for (int j = 0; j < cur_size; j++) {
      if (result[j]->getMatchEnvSize() < mng_ptr->max_env_num_per_vertex_) {
        VertexAPtr vertex (new VertexA(result[j]));
        if (vertex->addPreEnv(prev_match_envs[i], mng_ptr->max_env_num_per_peak_)) {
          vertex->trim();
          result.push_back(vertex);
        }
      }
    }
  }
  // add envelopes in current window 
  for (size_t i = 0; i < cur_match_envs.size(); i++) {
    cur_size = result.size();
    for (int j = 0; j < cur_size; j++) {
      if (result[j]->getMatchEnvSize() < mng_ptr->max_env_num_per_vertex_) {
        VertexAPtr vertex (new VertexA(result[j]));
        if (vertex->addCurEnv(cur_match_envs[i], mng_ptr->max_env_num_per_peak_)) {
          vertex->trim();
          result.push_back(vertex);
        }
      }
    }
  }
  return result;
}

    
// this is a copy of the previous two methods except that env_num is added
void VertexBase::addEmptyVertexB(FeatureMngPtr mng_ptr, VertexBPtrVec &result,
                                 DeconvDataPtr data, int cur_win, int env_num) {
  int prev_win = cur_win - 1;
  int bgn_peak = getBgnPeak(prev_win, data);
  int prev_win_peak_num = getWinPkNum(prev_win, data);
  int cur_win_peak_num = getWinPkNum(cur_win, data);
  VertexBPtr vertex(new VertexB(mng_ptr, bgn_peak, prev_win_peak_num, cur_win_peak_num, env_num));
  vertex->trim();
  result.push_back(vertex);
}

VertexBPtrVec VertexBase::getVertexBList(DeconvDataPtr data, int cur_win,
                                         MatchEnvPtrVec &prev_match_envs, 
                                         MatchEnvPtrVec &cur_match_envs, FeatureMngPtr mng_ptr) {
  VertexBPtrVec result;
  // add an empty vertex 
  addEmptyVertexB(mng_ptr, result, data, cur_win, mng_ptr->dp_env_num_);
  int cur_size = 1;
  // add envelopes in previous window 
  for (size_t i = 0; i < prev_match_envs.size(); i++) {
    cur_size = result.size();
    for (int j = 0; j < cur_size; j++) {
      if (result[j]->getMatchEnvSize() < mng_ptr->max_env_num_per_vertex_) {
        VertexBPtr vertex(new VertexB(result[j]));
        if (vertex->addPreEnv(prev_match_envs[i], mng_ptr->max_env_num_per_peak_)) {
          vertex->trim();
          result.push_back(vertex);
        }
      }
    }
  }
  // add envelopes in current window 
  for (size_t i = 0; i < cur_match_envs.size(); i++) {
    cur_size = result.size();
    for (int j = 0; j < cur_size; j++) {
      if (result[j]->getMatchEnvSize() < mng_ptr->max_env_num_per_vertex_) {
        VertexBPtr vertex(new VertexB(result[j]));
        if (vertex->addCurEnv(cur_match_envs[i], mng_ptr->max_env_num_per_peak_)) {
          vertex->trim();
          result.push_back(vertex);
        }
      }
    }
  }
  return result;
}

}
