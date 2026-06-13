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

#ifndef TOPPIC_TOPFD_DP_VERTEX_UTIL_HPP_
#define TOPPIC_TOPFD_DP_VERTEX_UTIL_HPP_

#include "topfd/deconv/deconv_data.hpp"
#include "topfd/dp/vertex_a.hpp"
#include "topfd/dp/vertex_b.hpp"

namespace toppic {

namespace vertex_util {

int getBgnPeak(int pre_win, const DeconvDataPtr& data);

int getWinPkNum(int win, const DeconvDataPtr& data);

void addEmptyVertexA(const DpParaPtr& dp_para_ptr, VertexAPtrVec& result,
                     const DeconvDataPtr& data, int cur_win);

VertexAPtrVec getVertexAList(const DeconvDataPtr& data, int cur_win,
                             const MatchEnvPtrVec& prev_match_envs,
                             const MatchEnvPtrVec& cur_match_envs,
                             const DpParaPtr& dp_para_ptr);

void addEmptyVertexB(const DpParaPtr& dp_para_ptr, VertexBPtrVec& result,
                     const DeconvDataPtr& data, int cur_win, int env_num);

VertexBPtrVec getVertexBList(const DeconvDataPtr& data, int cur_win,
                             MatchEnvPtrVec& prev_match_envs,
                             MatchEnvPtrVec& cur_match_envs,
                             const DpParaPtr& dp_para_ptr);
};  // namespace vertex_util

}  // namespace toppic
#endif
