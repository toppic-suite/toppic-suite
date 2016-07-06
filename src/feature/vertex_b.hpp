#ifndef PROT_FEATURE_VERTEX_B_HPP_
#define PROT_FEATURE_VERTEX_B_HPP_

#include "feature/vertex.hpp"

namespace prot {

class VertexB;
typedef std::shared_ptr<VertexB> VertexBPtr;

class VertexB : public Vertex {
 public:
  VertexB(FeatureMngPtr mng_ptr, int bgn_peak, 
          int pre_win_peak_num, int cur_win_peak_num, int env_num);
  VertexB(VertexBPtr ptr);

  bool addPreEnv(MatchEnvPtr env, int max_overlap);

  double getScoreB(int i) {return scores_[i];}

  int getPrevB(int i) {return prevs_[i];}

  void setScoreB(int i, double s) {scores_[i] = s;}

  void setPrevB(int i, int p) {prevs_[i] = p;}

 private:
    // dp scores for different number of envelopes 
  std::vector<double> scores_;
    // previous vertex id 
  std::vector<int> prevs_;
};

typedef std::vector<VertexBPtr> VertexBPtrVec;

}
#endif
