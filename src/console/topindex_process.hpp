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

#ifndef TOPPIC_TOPINDEX_PROCESS_HPP
#define TOPPIC_TOPINDEX_PROCESS_HPP

#include <string>
#include <map>
#include <vector>

#include "common/base/mod.hpp"
#include "seq/proteoform.hpp"

#include "filter/zeroptm/zero_ptm_filter_mng.hpp"
#include "filter/oneptm/one_ptm_filter_mng.hpp"
#include "filter/diag/diag_filter_mng.hpp"

namespace toppic {
    void TopIndexProcess(std::map<std::string, std::string> & arguments);
    void process(ZeroPtmFilterMngPtr zero_ptr, int thread_num);
    void process(OnePtmFilterMngPtr one_ptm_filter_mng_ptr, int thread_num);
    void process(DiagFilterMngPtr diag_filter_mng_ptr, int thread_num);
    
    void createIndexFiles(ProteoformPtrVec &raw_forms, int block_idx, ZeroPtmFilterMngPtr mng_ptr);
}  // namespace toppic

#endif
