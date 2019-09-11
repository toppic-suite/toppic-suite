//Copyright (c) 2014 - 2019, The Trustees of Indiana University.
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


#ifndef TOPPIC_DECONV_DP_VERTEX_HPP_
#define TOPPIC_DECONV_DP_VERTEX_HPP_

#include "deconv/dp/dp_para.hpp"
#include "deconv/env/match_env.hpp"
#include "deconv/env/env_peak_pair.hpp"

namespace toppic {

class Vertex;
typedef std::shared_ptr<Vertex> VertexPtr;

class Vertex {
 public:
  Vertex(DpParaPtr dp_para_ptr, int bgn_peak, int pre_win_peak_num, int cur_win_peak_num);

  Vertex(VertexPtr ptr);

  Vertex() {};

  virtual bool addPreEnv(MatchEnvPtr env, int max_overlap);
  bool addCurEnv(MatchEnvPtr env, int max_overlap);
  void trim();
  static bool checkConsist(VertexPtr pre, VertexPtr cur, int max_env_per_peak);

  int getMatchEnvSize() {return prec_match_envs_.size() + cur_match_envs_.size();}

  std::vector<int> getPeakUseCnts() {return peak_use_cnts_;}

  MatchEnvPtrVec getPreMatchEnvs() {return prec_match_envs_;}

  MatchEnvPtrVec getCurMatchEnvs() {return cur_match_envs_;}

  int getPreMatchEnvSize() {return prec_match_envs_.size();}

  static double getShareScr(VertexPtr pre, VertexPtr cur, double error_tole);

 protected:
  DpParaPtr dp_para_ptr_;
  int bgn_peak_;
  int peak_num_;
  int cur_bgn_pos_;

  std::vector<int> peak_use_cnts_;

  // for each envelope, there is a peak list 
  MatchEnvPtrVec prec_match_envs_;
  MatchEnvPtrVec cur_match_envs_;

  // There is a list of EnvPeakPairs for each peak 
  EnvPeakPairPtr2D pre_env_peak_pairs_;
  EnvPeakPairPtr2D cur_env_peak_pairs_;

  bool passDblIncrCheck(MatchEnvPtr env);

  bool passDblIncrCheck(MatchEnvPtrVec &env_list); 
};

typedef std::vector<VertexPtr> VertexPtrVec;

}
#endif
