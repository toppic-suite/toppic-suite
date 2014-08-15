#ifndef PROT_PTM_MNG_HPP_
#define PROT_PTM_MNG_HPP_

#include "prsm/prsm_para.hpp"

namespace prot {

class PtmMng {
 public :
  PtmMng(PrsmParaPtr prsm_para_ptr, int n_report, int shift_num,
         double align_max_shift, const std::string &input_file_ext, 
         const std::string &output_file_ext) {
    prsm_para_ptr_ = prsm_para_ptr;
    n_report_ = n_report;
    n_unknown_shift_ = shift_num;
    align_max_shift_ = align_max_shift;
    input_file_ext_ = input_file_ext;
    output_file_ext_ = output_file_ext;
  }
  
  PrsmParaPtr prsm_para_ptr_;

  /* parameters for ptm search */
  int n_report_ = 1;
  int n_unknown_shift_ =2;

  std::string input_file_ext_;
  std::string output_file_ext_;

  /* parameters for compute shift low memory */ 
  int ptm_fast_filter_scale_ = 100;
  int n_top_diagonals_ = 20;
  double min_double_gap=0.25;
  int min_diagonal_gap_ = (int)(ptm_fast_filter_scale_ * min_double_gap);

  /* parameters for diagonal generation */
  double extend_trunc_error_tolerance_ = 0.5;
  double align_prefix_suffix_shift_thresh_ = 300;

  /* parameters for ps alignment */
  double align_max_shift_ = 1000000;
  // remove min shift so that large errors in precursor mass can be treated 
  // a as small shift
  //double align_min_shift_ = 0.5;
  
  // shift thresh for penalty
  double align_large_shift_thresh_ = 300;
  // set panelty to 3
  double align_large_shift_panelty_ = 3;

  double refine_prec_step_width_ = 0.005;

};

typedef std::shared_ptr<PtmMng> PtmMngPtr;

} /* namespace prot */

#endif /* PTM_MNG_HPP_ */
