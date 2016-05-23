#ifndef PROT_DIAG_FILTER_MNG_HPP_
#define PROT_DIAG_FILTER_MNG_HPP_

#include "prsm/prsm_para.hpp"

namespace prot {

class DiagFilterMng {
 public:

  DiagFilterMng(PrsmParaPtr prsm_para_ptr, int filtering_result_num,
                int thread_num,
                const std::string &output_file_ext) {
    prsm_para_ptr_ = prsm_para_ptr;
    filter_result_num_ = filtering_result_num;
    thread_num_ = thread_num;
    output_file_ext_ = output_file_ext;
  }

  PrsmParaPtr prsm_para_ptr_;

  /** parameters for fast filteration */
  int max_proteoform_mass_ = 20000;

  //Candidate protein number for each spectrum
  size_t filter_result_num_ = 20;
  int db_block_size_ = 5000000;
  int filter_scale_ = 100;
  int thread_num_ = 1;

  std::string output_file_ext_;

};

typedef std::shared_ptr<DiagFilterMng> DiagFilterMngPtr;

} /* namespace tools */

#endif /* DIAG_FILTER_MNG_HPP_ */
