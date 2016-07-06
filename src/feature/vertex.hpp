#ifndef PROT_FEATURE_VERTEX_HPP_
#define PROT_FEATURE_VERTEX_HPP_

#include "feature/feature_mng.hpp"
#include "feature/match_env.hpp"
#include "feature/env_peak_pair.hpp"

namespace prot {

class Vertex;
typedef std::shared_ptr<Vertex> VertexPtr;

class Vertex {
 public:
  Vertex(FeatureMngPtr mng_ptr, int bgn_peak, int pre_win_peak_num, int cur_win_peak_num);

  Vertex(VertexPtr ptr);

  bool addPreEnv(MatchEnvPtr env, int max_overlap);
  bool addCurEnv(MatchEnvPtr env, int max_overlap);
  void trim();
  bool checkConsist(VertexPtr pre, VertexPtr cur, int max_env_per_peak);

  int getMatchEnvSize() {return prec_match_envs_.size() + cur_match_envs_.size();}

  std::vector<int> getPeakUseCnts() {return peak_use_cnts_;}

  MatchEnvPtrVec getPreMatchEnvs() {return prec_match_envs_;}

  MatchEnvPtrVec getCurMatchEnvs() {return cur_match_envs_;}

  int getPreMatchEnvSize() {return prec_match_envs_.size();}

  static double getShareScr(VertexPtr pre, VertexPtr cur, double error_tole);

 private:
  FeatureMngPtr mng_ptr_;
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
