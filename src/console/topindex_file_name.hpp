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

//#ifndef TOPPIC_CONSOLE_TOPINDEX_FILE_NAME_HPP_
//#define TOPPIC_CONSOLE_TOPINDEX_FILE_NAME_HPP_

#include <vector>
#include <map>
#include <vector>

#include "prsm/prsm_para.hpp"

namespace toppic{
class TopIndexFileName{
    public:
    TopIndexFileName(){};
    std::string gene_file_name(PrsmParaPtr prsm_para_ptr);

    std::vector<std::string> zero_ptm_file_vec{"zero_ptm_term_index_", "zero_ptm_diag_index_", 
    "zero_ptm_rev_term_index_", "zero_ptm_rev_diag_index_"};//file name vector

    std::vector<std::string> one_ptm_file_vec{"one_ptm_term_index_", "one_ptm_diag_index_", 
    "one_ptm_rev_term_index_", "one_ptm_rev_diag_index_"};//file name vector

    std::vector<std::string> multi_ptm_file_vec{"multi_ptm_index_"};//file name vector

    std::map<std::string, std::string> prot_mod_map = {
        {"NONE", "N"}, {"NME", "NME"}, {"NME_ACETYLATION", "NMEA"}, {"M_ACETYLATION", "MA"}
    };
    std::map<std::string, std::string> search_type_map = {
        {"TARGET", "no_decoy"}, {"TARGET+DECOY", "decoy"}};
    };
}
