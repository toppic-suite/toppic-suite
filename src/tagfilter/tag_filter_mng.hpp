
#ifndef PROT_TAG_FILTER_MNG
#define PROT_TAG_FILTER_MNG

#include <string>

#include "prsm/prsm_para.hpp"

namespace prot {

class TagFilterMng {
 public:

  TagFilterMng(PrsmParaPtr prsm_para_ptr, int filtering_result_num,
               const std::string &output_file_ext) {
    prsm_para_ptr_ = prsm_para_ptr;
    filter_result_num_ = filtering_result_num;
    output_file_ext_ = output_file_ext;
  }

  PrsmParaPtr prsm_para_ptr_;
  
  int filter_result_num_ = 20;

  int prec_error_ = 2;

  std::string output_file_ext_;

};

typedef std::shared_ptr<TagFilterMng> TagFilterMngPtr;

}
#endif
