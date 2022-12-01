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

#ifndef TOPPIC_FILTER_MNG_ONE_PTM_FILTER_MNG_HPP_
#define TOPPIC_FILTER_MNG_ONE_PTM_FILTER_MNG_HPP_

#include <boost/thread.hpp>

#include "para/prsm_para.hpp"

namespace toppic {

class OnePtmFilterMng {
 public:
  // for index generation
  OnePtmFilterMng(PrsmParaPtr prsm_para_ptr,
                  const std::string & index_file_para,
                  const std::string & output_file_ext,
                  int thread_num);

  // for one shift filter in TopPIC
  OnePtmFilterMng(PrsmParaPtr prsm_para_ptr,
                  const std::string & index_file_para,
                  const std::string & output_file_ext,
                  int thread_num,
                  double min_shift,
                  double max_shift); 

  // for One shift filter in TopMG
  OnePtmFilterMng(PrsmParaPtr prsm_para_ptr,
                  const std::string & index_file_para,
                  const std::string & output_file_ext,
                  int thread_num,
                  double min_shift,
                  double max_shift,  
                  const std::string & residueModFileName, 
                  int var_num);

  std::string getIndexFilePara() {return index_file_para_;}

  PrsmParaPtr prsm_para_ptr_;


  // parameters for fast filteration
  int max_proteoform_mass_ = 50000;

  // Candidate protein number for each spectrum
  unsigned int comp_num_ = 10;
  unsigned int pref_suff_num_ = 5;
  unsigned int inte_num_ = 10;
  
  // reported n term or c term shift number
  unsigned int comp_cand_shift_num_ = 1;
  unsigned int cand_shift_num_ = 10;
  int filter_scale_ = 100;

  int threshold_ = 4;

  std::string index_file_para_;

  std::string output_file_ext_;

  int thread_num_;

  double min_shift_= -500;
  double max_shift_ = 500;

  std::string residueModFileName_ = "";

  int var_num_ = 0;

  int n_spec_block_ = 0;

  boost::mutex mutex_;

  std::vector<int> cnts_;

  //file name vector
  std::vector<std::string> one_ptm_file_vec_{"one_ptm_term_index", "one_ptm_diag_index", 
    "one_ptm_rev_term_index", "one_ptm_rev_diag_index"};
};

typedef std::shared_ptr<OnePtmFilterMng> OnePtmFilterMngPtr;

}  // namespace toppic

#endif /* ONE_PTM_FILTER_MNG_HPP_ */
