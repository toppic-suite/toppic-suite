// Copyright (c) 2014 - 2016, The Trustees of Indiana University.
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
