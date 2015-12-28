#ifndef PROT_PTM_SEARCH_MNG_HPP_
#define PROT_PTM_SEARCH_MNG_HPP_

#include "prsm/prsm_para.hpp"
#include "oneptmsearch/ps_align_para.hpp"

namespace prot {

class PtmSearchMng {
 public :
  PtmSearchMng(PrsmParaPtr prsm_para_ptr, int n_report, int shift_num,
               double align_max_shift, const std::string &input_file_ext, 
               const std::string &output_file_ext) {
    prsm_para_ptr_ = prsm_para_ptr;
    n_report_ = n_report;
    input_file_ext_ = input_file_ext;
    output_file_ext_ = output_file_ext;
    align_para_ptr_ = PsAlignParaPtr(new PsAlignPara(shift_num, align_max_shift));
  }
  
  int getShiftNum() {return align_mng_ptr_->n_unknown_shift_;}
  
  PrsmParaPtr prsm_para_ptr_;

  /* parameters for ptm search */
  int n_report_ = 1;

  std::string input_file_ext_;
  std::string output_file_ext_;

  /* parameters for compute shift low memory */ 
  int ptm_fast_filter_scale_ = 100;
  int n_top_diagonals_ = 20;
  double min_double_gap = 0.25;
  int min_diagonal_gap_ = (int)(ptm_fast_filter_scale_ * min_double_gap);

  /* parameters for diagonal generation */
  double extend_trunc_error_tolerance_ = 0.5;
  double align_prefix_suffix_shift_thresh_ = 300;

  PsAlignMngPtr align_para_ptr_;

};

typedef std::shared_ptr<PtmMng> PtmMngPtr;

} /* namespace prot */

#endif /* PTM_SEARCH_MNG_HPP_ */
