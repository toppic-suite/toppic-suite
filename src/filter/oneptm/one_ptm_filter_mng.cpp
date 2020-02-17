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

#include "filter/oneptm/one_ptm_filter_mng.hpp"

namespace toppic {

OnePtmFilterMng::OnePtmFilterMng(PrsmParaPtr prsm_para_ptr,
                                 const std::string & output_file_ext,
                                 int thread_num,
                                 const std::string & residueModFileName,
                                 int var_num):
    prsm_para_ptr_(prsm_para_ptr),
    output_file_ext_(output_file_ext),
    thread_num_(thread_num),
    residueModFileName_(residueModFileName),
    var_num_(var_num) {}

}  // namespace toppic

