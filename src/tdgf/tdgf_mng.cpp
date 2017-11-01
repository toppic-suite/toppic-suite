//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#include "tdgf/tdgf_mng.hpp"

namespace prot {

TdgfMng::TdgfMng(PrsmParaPtr prsm_para_ptr, 
                 int shift_num, double max_ptm_mass, bool use_gf,
                 bool variable_ptm, int thread_num,
                 const std::string &input_file_ext, 
                 const std::string &output_file_ext):
    input_file_ext_(input_file_ext),
    output_file_ext_(output_file_ext),
    prsm_para_ptr_(prsm_para_ptr),
    use_gf_(use_gf),
    variable_ptm_(variable_ptm),
    max_ptm_mass_(max_ptm_mass),
    unexpected_shift_num_(shift_num),
    thread_num_(thread_num) {
    }

}
