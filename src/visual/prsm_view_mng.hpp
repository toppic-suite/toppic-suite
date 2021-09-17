//Copyright (c) 2014 - 2020, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.


#ifndef TOPPIC_VISUAL_PRSM_VIEW_MNG_HPP_
#define TOPPIC_VISUAL_PRSM_VIEW_MNG_HPP_

#include <string>

#include "common/util/logger.hpp"
#include "common/util/file_util.hpp"
#include "para/prsm_para.hpp"

namespace toppic {

class PrsmViewMng {
 public:
  PrsmViewMng(PrsmParaPtr prsm_para_ptr,
              const std::string & resource_dir,
              const std::string & fname_suffix):
      prsm_para_ptr_(prsm_para_ptr) {
        std::string spectrum_file_name = prsm_para_ptr_->getSpectrumFileName();
        std::string base_name = file_util::basename(spectrum_file_name);
        xml_path_ = base_name + "_" + fname_suffix + "_xml";
        html_path_ = base_name.substr(0, base_name.length() - 4) + "_html" 
            + file_util::getFileSeparator() + fname_suffix;
        LOG_DEBUG("html path " << html_path_);
        resource_dir_ = resource_dir;
        min_mass_ = prsm_para_ptr_->getSpParaPtr()->getMinMass();
      }

  PrsmParaPtr prsm_para_ptr_;

  std::string html_path_;

  std::string xml_path_;

  std::string resource_dir_;

  int decimal_point_num_ = 2;

  int precise_point_num_ = 4;

  double min_mass_;
};

typedef std::shared_ptr<PrsmViewMng> PrsmViewMngPtr;

} /* namespace toppic */


#endif /* TOPPIC_VISUAL_PRSM_VIEW_MNG_HPP_ */
