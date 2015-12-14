#ifndef PROT_TDGF_MNG_HPP_
#define PROT_TDGF_MNG_HPP_

#include "prsm/prsm_para.hpp"

namespace prot {

class TdgfMng {
 public:
  TdgfMng(PrsmParaPtr prsm_para_ptr, int shift_num, double max_ptm_mass, bool use_gf, 
          bool variable_ptm, const std::string &input_file_ext, 
          const std::string &output_file_ext);

  /*********************************
   * Common parameters
   * *******************************/
  std::string input_file_ext_;
  std::string output_file_ext_;

  PrsmParaPtr prsm_para_ptr_;

  /** Prsm filter */
  int comp_evalue_min_match_frag_num_ = 4;

  bool use_gf_ = false;

  bool variable_ptm_ = false;

  /**********************************
   * Tdgf Table parameters 
   * ********************************/


  /**********************************
   * Tdgf parameters 
   * ********************************/
  /* max ptm mass is used in the function for counting sequence numbers*/
  double max_ptm_mass_ = 1000000;

  /** do tdgf computation if poisson report evalue > 10^-8 
   * or match frag num < 25 */
  double computation_evalue_cutoff = 0.00000001;
  int computation_frag_num_cutoff = 25;

  /** dp table */
  // number of mass shift
  int unexpected_shift_num_ = 2;
  double convert_ratio_ = 274.335215;
  double max_prec_mass_ = 51000;
  int max_table_height_ = 128;
  int min_height_ = 10;

};

typedef std::shared_ptr<TdgfMng> TdgfMngPtr;

}
#endif
