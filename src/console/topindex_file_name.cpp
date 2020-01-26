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

#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>  

#include "prsm/prsm_para.hpp"
#include "console/topindex_file_name.hpp"

namespace toppic{

std::string TopIndexFileName::gene_file_name(PrsmParaPtr prsm_para_ptr){
    std::string fixed_mod = prsm_para_ptr->getFixedMod();
    std::string error_tol = prsm_para_ptr->getErrorTolerance();
    std::string decoy = search_type_map.find(prsm_para_ptr->getSearchType())->second;
    std::string activation = prsm_para_ptr->getActivation();
    std::string prot_mod = prsm_para_ptr->getProtMod();

    std::vector<std::string> prot_mod_vec;
    std::stringstream sstream(prot_mod);
    std::string final_prot_mod;

    while(sstream.good()){
        std::string substring;
        getline(sstream, substring, ',');
        prot_mod_vec.push_back(substring);
    }
    for (size_t i = 0; i < prot_mod_vec.size(); i++){
        final_prot_mod = final_prot_mod + "_" + prot_mod_map.find(prot_mod_vec[i])->second;
    }
    std::string paraInfo =  fixed_mod + "_" + final_prot_mod + "_" + activation + "_" + error_tol + "_" + decoy;
    
    return paraInfo;
    }
}