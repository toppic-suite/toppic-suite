//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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


#ifndef PROT_ZERO_PTM_SEARCH_MNG_HPP_
#define PROT_ZERO_PTM_SEARCH_MNG_HPP_

#include <string>

#include "prsm/prsm_para.hpp"

namespace prot {

class ZeroPtmSearchMng {
 public:
  ZeroPtmSearchMng(PrsmParaPtr prsm_para_ptr,
                   const std::string & input_file_ext,
                   const std::string & output_file_ext):
      prsm_para_ptr_(prsm_para_ptr),
      input_file_ext_(input_file_ext),
      output_file_ext_(output_file_ext) {}

  PrsmParaPtr prsm_para_ptr_;

  /** parameters for zero ptm search */

  /** zero ptm fast filtering */
  int zero_ptm_filter_result_num_ = 20;
  /** number of reported Prsms for each spectrum */
  int report_num_ = 1;

  /** recalibration is used in ZeroPtmSlowMatch */
  bool   do_recalibration_ = false;
  double recal_ppo_ = 0.000015;  // 15 ppm
  bool   ms_one_ms_two_same_recal_ = true;

  std::string input_file_ext_;

  std::string output_file_ext_;
};

typedef std::shared_ptr<ZeroPtmSearchMng> ZeroPtmSearchMngPtr;

}  // namespace prot

#endif
