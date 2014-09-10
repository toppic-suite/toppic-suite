#ifndef PROT_PRSM_VIEW_MNG_HPP_
#define PROT_PRSM_VIEW_MNG_HPP_

#include "base/file_util.hpp"
#include "prsm/prsm_para.hpp"

namespace prot {

class PrsmViewMng {
 public:
  PrsmViewMng(PrsmParaPtr prsm_para_ptr, std::string exec_dir);

  PrsmParaPtr prsm_para_ptr_;

  std::string html_path_;
  std::string xml_path_;
  std::string executive_dir_;

  int decimal_point_num_ = 2;
  int precise_point_num_ = 4;
  double min_mass_;
};

typedef std::shared_ptr<PrsmViewMng> PrsmViewMngPtr;

inline PrsmViewMng::PrsmViewMng(PrsmParaPtr prsm_para_ptr, std::string exec_dir) {
  prsm_para_ptr_ = prsm_para_ptr;
  std::string spectrum_file_name = prsm_para_ptr_->getSpectrumFileName();

  xml_path_ = basename(spectrum_file_name) + "_xml";
  html_path_ = basename(spectrum_file_name) + "_html";
  executive_dir_ = exec_dir;
  min_mass_ = prsm_para_ptr_->getSpParaPtr()->getMinMass();
}

} /* namespace prot */


#endif /* PROT_PRSM_VIEW_MNG_HPP_ */
