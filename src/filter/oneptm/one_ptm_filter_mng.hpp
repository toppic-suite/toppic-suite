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

#ifndef TOPPIC_FILTER_ONE_PTM_ONE_PTM_FILTER_MNG_HPP_
#define TOPPIC_FILTER_ONE_PTM_ONE_PTM_FILTER_MNG_HPP_

#include <atomic>

#include "prsm/prsm_para.hpp"

namespace toppic {

class OnePtmFilterMng {
 public:
  OnePtmFilterMng(PrsmParaPtr prsm_para_ptr,
                  const std::string & output_file_ext,
                  int thread_num,
                  const std::string & residueModFileName = "",
                  int var_num = 0);

  PrsmParaPtr prsm_para_ptr_;

  // parameters for fast filteration
  int max_proteoform_mass_ = 50000;

  // Candidate protein number for each spectrum
  unsigned int comp_num_ = 10;
  unsigned int pref_suff_num_ = 5;
  unsigned int inte_num_ = 10;
  unsigned int shift_num_ = 10;
  int filter_scale_ = 100;

  std::string output_file_ext_;

  int thread_num_;

  std::string residueModFileName_;

  int var_num_;

  std::atomic<int> cnt_;
  int n_spec_block_ = 0;
};

typedef std::shared_ptr<OnePtmFilterMng> OnePtmFilterMngPtr;

}  // namespace toppic

#endif /* ONE_PTM_FILTER_MNG_HPP_ */
