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
