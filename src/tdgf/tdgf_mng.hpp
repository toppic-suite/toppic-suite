#ifndef PROT_TDGF_MNG_HPP_
#define PROT_TDGF_MNG_HPP_

#include "prsm/prsm_para.hpp"

namespace prot {

class TdgfMng {
 public:
  TdgfMng(PrsmParaPtr prsm_para_ptr, int shift_num, double max_shift_value, 
          std::string input_file_ext, std::string output_file_ext) {
    prsm_para_ptr_ = prsm_para_ptr,
    max_shift_value_ = max_shift_value;
    unexpected_shift_num_ = shift_num;
    input_file_ext_ = input_file_ext;
    output_file_ext_ = output_file_ext;
  }

  PrsmParaPtr prsm_para_ptr_;

  /** parameters for tdmg */
  double max_shift_value_ = 1000000;

  /** Prsm filter */
  int comp_evalue_min_match_frag_num_ = 4;

  /** dp table */
  // number of mass shift
  int unexpected_shift_num_ = 2;
  double double_to_int_constant_ = 274.335215;
  double max_sp_prec_mass_ = 51000;
  int max_table_height_ = 128;
  int min_height_ = 10;

  /* Semi alignment type determination */
  /* a prsm with a shift < 300 at n-terminus is treated as a prefix */
  double prefix_suffix_shift_thresh_ = 300;

  /**
   * Postprocessing: adjustment makes it more conservative to identify Prsms
   * with multiple shifts
   */
  // the following three values should be adjusted
  double multiple_shift_adjustment_ = 4;
  double multiple_shfit_adjustment_base_ = 10;
  double min_adjustment_ = 100;

  std::string input_file_ext_;
  std::string output_file_ext_;
};

typedef std::shared_ptr<TdgfMng> TdgfMngPtr;

}
#endif
