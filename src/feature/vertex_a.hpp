#ifndef PROT_FEATURE_VERTEX_A_HPP_
#define PROT_FEATURE_VERTEX_A_HPP_

#include "feature/vertex.hpp"

namespace prot {

class VertexA;
typedef std::shared_ptr<VertexA> VertexAPtr;

class VertexA : public Vertex {
 public:
  VertexA(FeatureMngPtr mng_ptr, int bgn_peak, 
          int pre_win_peak_num, int cur_win_peak_num);
  VertexA(VertexAPtr ptr);

  bool addPreEnv(MatchEnvPtr env, int max_overlap);

  double getThisScr() {return this_score_;}

  double getScrA() {return score_;}

  void setScrA(double score) {score_ = score;}

  int getPreA() {return prev_vertex_;}

  void setPreA(int prev) {prev_vertex_ = prev;}

 private:
  // the sum of score_s of envelopes in previous window 
  double this_score_;
  // current score_ for dp, only the score_s of pre_match_env are included. The
  // reason is that we can only determine the score_ for sharing model when
  // envelopes in the next window are determined.
  double score_;
  // previous vertex for backtracking 
  int prev_vertex_;
};

typedef std::vector<VertexAPtr> VertexAPtrVec;
typedef std::vector<VertexAPtrVec> VertexAPtr2D;

}
#endif
