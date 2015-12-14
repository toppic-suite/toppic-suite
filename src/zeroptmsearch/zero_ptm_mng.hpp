#ifndef PROT_ZERO_PTM_MNG_HPP_
#define PROT_ZERO_PTM_MNG_HPP_

#include <string>

#include "prsm/prsm_para.hpp"

namespace prot {

class ZeroPtmMng {
 public:
  ZeroPtmMng(PrsmParaPtr prsm_para_ptr, 
             const std::string & input_file_ext,
             const std::string & output_file_ext): 
      prsm_para_ptr_(prsm_para_ptr), 
      input_file_ext_(input_file_ext),
      output_file_ext_(output_file_ext) {
      }

  PrsmParaPtr prsm_para_ptr_;

  /** parameters for zero ptm search */

  /** zero ptm fast filtering */
  int zero_ptm_filter_result_num_ = 20;
  /** number of reported Prsms for each spectrum */
  int report_num_ = 1;

  /** recalibration is used in ZeroPtmSlowMatch */
  bool   do_recalibration_ = false;
  double recal_ppo_ = 0.000015; // 15 ppm
  bool   ms_one_ms_two_same_recal_ = true;

  std::string input_file_ext_;

  std::string output_file_ext_;
};

typedef std::shared_ptr<ZeroPtmMng> ZeroPtmMngPtr;

} /* namespace_prot */

#endif
