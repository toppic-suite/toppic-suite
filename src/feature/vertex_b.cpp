#include "feature/vertex_b.hpp"

namespace prot {

VertexB::VertexB(FeatureMngPtr mng_ptr, int bgn_peak, 
                 int pre_win_peak_num, int cur_win_peak_num, int env_num) {
  // zero envelope is considered, so we have env_num + 1 
  mng_ptr_ = mng_ptr;
  bgn_peak_ = bgn_peak;
  peak_num_ = pre_win_peak_num + cur_win_peak_num;
  cur_bgn_pos_ = pre_win_peak_num;
  peak_use_cnts_.resize(peak_num_, 0);

  pre_env_peak_pairs_.resize(peak_num_); 
  cur_env_peak_pairs_.resize(peak_num_); 
  scores_.resize(env_num + 1, 0);
  prevs_.resize(env_num + 1, -1);
}

VertexB::VertexB(VertexBPtr ptr) {
  mng_ptr_= ptr->mng_ptr_;
  bgn_peak_ = ptr->bgn_peak_;
  peak_num_ = ptr->peak_num_;
  cur_bgn_pos_ = ptr->cur_bgn_pos_;

  peak_use_cnts_ = ptr->peak_use_cnts_;

  prec_match_envs_ = ptr->prec_match_envs_;
  cur_match_envs_ = ptr->cur_match_envs_;

  pre_env_peak_pairs_ = ptr->pre_env_peak_pairs_;
  cur_env_peak_pairs_ = ptr->cur_env_peak_pairs_;
  scores_ = ptr->scores_;
  prevs_ = ptr->prevs_;
}

}
