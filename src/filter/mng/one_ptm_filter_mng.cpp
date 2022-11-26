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

#include "filter/mng/one_ptm_filter_mng.hpp"

namespace toppic {

OnePtmFilterMng::OnePtmFilterMng(PrsmParaPtr prsm_para_ptr,
                                 const std::string & index_file_para,
                                 const std::string & output_file_ext,
                                 int thread_num): 
    prsm_para_ptr_(prsm_para_ptr),
    index_file_para_(index_file_para),
    output_file_ext_(output_file_ext),
    thread_num_(thread_num) {}

  // for one shift filter in TopPIC
OnePtmFilterMng::OnePtmFilterMng(PrsmParaPtr prsm_para_ptr,
                                 const std::string & index_file_para,
                                 const std::string & output_file_ext,
                                 int thread_num,
                                 double min_shift,
                                 double max_shift):  
  prsm_para_ptr_(prsm_para_ptr),
  index_file_para_(index_file_para),
  output_file_ext_(output_file_ext),
  thread_num_(thread_num),
  min_shift_(min_shift),
  max_shift_(max_shift) {}

  // for One shift filter in TopMG
OnePtmFilterMng::OnePtmFilterMng(PrsmParaPtr prsm_para_ptr,
                                 const std::string & index_file_para,
                                 const std::string & output_file_ext,
                                 int thread_num,
                                 double min_shift,
                                 double max_shift,  
                                 const std::string & residueModFileName, 
                                 int var_num):
    prsm_para_ptr_(prsm_para_ptr),
    index_file_para_(index_file_para),
    output_file_ext_(output_file_ext),
    thread_num_(thread_num),
    min_shift_(min_shift),
    max_shift_(max_shift),
    residueModFileName_(residueModFileName),
    var_num_(var_num) {}

}  // namespace toppic

