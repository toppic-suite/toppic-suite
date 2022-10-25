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

#ifndef TOPPIC_FILTER_MNG_ZERO_PTM_FILTER_MNG_HPP_
#define TOPPIC_FILTER_MNG_ZERO_PTM_FILTER_MNG_HPP_

#include <boost/thread.hpp>

#include "para/prsm_para.hpp"

namespace toppic {

class ZeroPtmFilterMng {
 public:
  ZeroPtmFilterMng(PrsmParaPtr prsm_para_ptr,
                   const std::string &index_file_para,
                   int thread_num,
                   const std::string &output_file_ext);

  PrsmParaPtr getPrsmPtr(){return prsm_para_ptr_;}

  std::string getIndexFilePara() {return index_file_para_;}

  PrsmParaPtr prsm_para_ptr_;

  std::string index_file_para_;

  //Parameters for fast filteration 
  int max_proteoform_mass_ = 100000;

  //Candidate protein number for each spectrum
  unsigned int comp_num_ = 5;
  unsigned int pref_suff_num_ = 5;
  unsigned int inte_num_ = 10;
  int filter_scale_ = 100;

  int thread_num_ = 1;

  int n_spec_block_ = 0;

  std::string output_file_ext_;

  boost::mutex mutex_;

  std::vector<int> cnts_;
  
  std::vector<std::string> zero_ptm_file_vec_{"zero_ptm_term_index", "zero_ptm_diag_index", 
    "zero_ptm_rev_term_index", "zero_ptm_rev_diag_index"};
};

typedef std::shared_ptr<ZeroPtmFilterMng> ZeroPtmFilterMngPtr;

} 

#endif 
