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

#ifndef TOPPIC_FILTER_DIAG_DIAG_FILTER_MNG_HPP_
#define TOPPIC_FILTER_DIAG_DIAG_FILTER_MNG_HPP_

#include "prsm/prsm_para.hpp"

namespace toppic {

class DiagFilterMng {
 public:
  DiagFilterMng(PrsmParaPtr prsm_para_ptr,
                int filtering_result_num,
                int thread_num,
                const std::string &output_file_ext,
                const std::string & residueModFileName = "",
                int var_num = 0);

  PrsmParaPtr prsm_para_ptr_;

  /** parameters for fast filteration */
  int max_proteoform_mass_ = 40000;

  // Candidate protein number for each spectrum
  size_t filter_result_num_ = 20;
  int db_block_size_ = 5000000;
  int filter_scale_ = 100;
  int thread_num_ = 1;

  std::string output_file_ext_;
  std::string residueModFileName_;

  int var_num_;

    std::vector<std::string> file_names{"toppic_multi_ptm_complete", "toppic_multi_ptm_prefix", "toppic_multi_ptm_suffix", "toppic_multi_ptm_internal"};
};

typedef std::shared_ptr<DiagFilterMng> DiagFilterMngPtr;

}  // namespace toppic

#endif /* DIAG_FILTER_MNG_HPP_ */
