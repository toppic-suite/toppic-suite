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

#include "filter/mng/zero_ptm_filter_mng.hpp"

namespace toppic {

ZeroPtmFilterMng::ZeroPtmFilterMng(PrsmParaPtr prsm_para_ptr,
                                   int thread_num,
                                   const std::string & output_file_ext):
    prsm_para_ptr_(prsm_para_ptr),
    thread_num_(thread_num),
    output_file_ext_(output_file_ext) {}
} 
