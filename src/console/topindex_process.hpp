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
#include "common/thread/simple_thread_pool.hpp"

#include "seq/proteoform.hpp"
#include "seq/db_block.hpp"

#include "filter/zeroptm/mass_zero_ptm_filter.hpp"
#include "filter/zeroptm/zero_ptm_filter_mng.hpp"
#include "filter/oneptm/one_ptm_filter_mng.hpp"
#include "filter/diag/diag_filter_mng.hpp"

#include "prsm/prsm_para.hpp"

namespace toppic {
    void TopIndexProcess(std::map<std::string, std::string> & arguments);
    void index_process();
    std::string gene_file_name(PrsmParaPtr prsm_para_ptr);

    std::vector<std::string> zero_ptm_file_vec{"zero_ptm_term_index", "zero_ptm_diag_index", 
    "zero_ptm_rev_term_index", "zero_ptm_rev_diag_index"};//file name vector

    std::map<std::string, std::string> prot_mod_map = {
        {"NONE", "N"}, {"NME", "NME"}, {"NME_ACETYLATION", "NMEA"}, {"M_ACETYLATION", "MA"}
    };
    std::map<std::string, std::string> search_type_map = {
        {"TARGET", ""}, {"TARGET+DECOY", "decoy"}};

}  // namespace toppic

#endif
