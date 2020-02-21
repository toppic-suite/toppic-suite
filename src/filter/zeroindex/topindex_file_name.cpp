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
#include <sstream> 

#include "prsm/prsm_para.hpp"
#include "filter/zeroindex/topindex_file_name.hpp"

namespace toppic{

std::string TopIndexFileName::geneFileName(PrsmParaPtr prsm_para_ptr){
    std::string fixed_mod = prsm_para_ptr->getFixedMod();
    std::string error_tol = prsm_para_ptr->getErrorTolerance();
    std::string decoy = search_type_map_.find(prsm_para_ptr->getSearchType())->second;
    std::string activation = prsm_para_ptr->getActivation();
    std::string prot_mod = prsm_para_ptr->getProtMod();

    std::vector<std::string> prot_mod_vec;
    std::stringstream sstream(prot_mod);
    std::string final_prot_mod = "";

    while(sstream.good()){
        std::string substring;
        getline(sstream, substring, ',');
        prot_mod_vec.push_back(substring);
    }
   
    for (size_t i = 0; i < prot_mod_vec.size(); i++){
        std::string mod = prot_mod_map_.find(prot_mod_vec[i])->second;
        if (final_prot_mod != ""){
           final_prot_mod = final_prot_mod + "_" + mod;     
        }else{
            final_prot_mod = final_prot_mod + mod;
        }
    }

    std::vector<std::string>para_vec;//to determine if "_" is needed in between
    std::string para_info;

    para_vec.push_back(fixed_mod);
    para_vec.push_back(final_prot_mod);
    para_vec.push_back(activation);
    para_vec.push_back(error_tol);
    para_vec.push_back(decoy);

    for (size_t t = 0; t < para_vec.size(); t++){
        if (para_vec[t] != ""){
            //if the parameter is used
            para_info = para_info + "_" + para_vec[t];
        }
    }
    return para_info;   
    }
}
