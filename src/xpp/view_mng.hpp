/*
 * view_mng.hpp
 *
 *  Created on: Mar 25, 2014
 *      Author: xunlikun
 */

#ifndef VIEW_MNG_HPP_
#define VIEW_MNG_HPP_

#include "base/file_util.hpp"
#include "prsm/prsm_para.hpp"

namespace prot {

class ViewMng {
 public:
  ViewMng(PrsmParaPtr prsm_para_ptr, std::string exec_dir) {
    prsm_para_ptr_ = prsm_para_ptr;
    std::string spectrum_file_name = prsm_para_ptr_->getSpectrumFileName();

    xml_path_ = basename(spectrum_file_name) + "_xml";
    html_path_ = basename(spectrum_file_name) + "_html";
    executive_dir_ = exec_dir;
  }

  PrsmParaPtr prsm_para_ptr_;

  std::string html_path_;
  std::string xml_path_;
  std::string executive_dir_;
};

typedef std::shared_ptr<ViewMng> ViewMngPtr;
} /* namespace prot */

#endif /* VIEW_MNG_HPP_ */
