/*
 * ptm_mng.hpp
 *
 *  Created on: Dec 19, 2013
 *      Author: xunlikun
 */

#ifndef PROT_PTM_MNG_HPP_
#define PROT_PTM_MNG_HPP_

#include "prsm/prsm_para.hpp"

namespace prot {

class PtmMng {
 public :
  PtmMng(PrsmParaPtr prsm_para_ptr, int n_report, int shift_num,
         std::string input_file_ext, std::string output_file_ext) {
    prsm_para_ptr_ = prsm_para_ptr;
    n_report_ = n_report;
    n_unknown_shift_ = shift_num;
    input_file_ext_ = input_file_ext;
    output_file_ext_ = output_file_ext;
  }
  
  PrsmParaPtr prsm_para_ptr_;

  /* parameters for ptm search */
  int n_report_ = 1;
  int n_unknown_shift_ =2;
  int n_known_shift_ = 0;

  int ptm_fast_filter_scale_ = 100;
  int n_top_diagonals_ = 20;
  double min_double_gap=0.25;
  int min_diagonal_gap_ = (int)(ptm_fast_filter_scale_ * min_double_gap);
  double extend_diagonal_error_tolerance_ = 0.5;
  double test_term_match_error_tolerance_ = 0.05;
  double align_prefix_suffix_shift_thresh_ = 300;

  double align_min_gap = 0.5;
  double large_shift_thresh = 300;
  double large_shift_panelty = 0;

  double adjust_prec_step_width_ = 0.005;

  std::string input_file_ext_;
  std::string output_file_ext_;
};

typedef std::shared_ptr<PtmMng> PtmMngPtr;

} /* namespace prot */

#endif /* PTM_MNG_HPP_ */
