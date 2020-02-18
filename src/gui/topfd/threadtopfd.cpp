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

#include "topfd/common/topfd_para.hpp"
#include "topfd/common/topfd_process.hpp"
#include "gui/topfd/threadtopfd.hpp"

void ThreadTopFD::setPar(toppic::TopfdParaPtr para_ptr, 
                         const std::vector<std::string> & spec_file_lst) {
  para_ptr_ = para_ptr;
  spec_file_lst_ = spec_file_lst;
}


void ThreadTopFD::run() {
  toppic::topfd_process::process(para_ptr_, spec_file_lst_);
}
