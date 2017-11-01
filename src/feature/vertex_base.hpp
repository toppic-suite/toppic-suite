//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#ifndef PROT_FEATURE_VERTEX_BASE_HPP_
#define PROT_FEATURE_VERTEX_BASE_HPP_

#include "feature/feature_mng.hpp"
#include "feature/deconv_data.hpp"
#include "feature/vertex_a.hpp"
#include "feature/vertex_b.hpp"

namespace prot {

class VertexBase {
 public:
  static int getBgnPeak(int pre_win, DeconvDataPtr data);

  static int getWinPkNum(int win, DeconvDataPtr data);

  static void addEmptyVertexA(FeatureMngPtr mng_ptr, VertexAPtrVec &result, 
                              DeconvDataPtr data, int cur_win);

  static VertexAPtrVec getVertexAList(DeconvDataPtr data, int cur_win, MatchEnvPtrVec prev_match_envs, 
                                      MatchEnvPtrVec cur_match_envs, FeatureMngPtr mng_ptr);

  static void addEmptyVertexB(FeatureMngPtr mng_ptr, VertexBPtrVec &result,
                              DeconvDataPtr data, int cur_win, int env_num);

  static VertexBPtrVec getVertexBList(DeconvDataPtr data, int cur_win,
                                      MatchEnvPtrVec &prev_match_envs, 
                                      MatchEnvPtrVec &cur_match_envs, FeatureMngPtr mng_ptr);
};

}
#endif
