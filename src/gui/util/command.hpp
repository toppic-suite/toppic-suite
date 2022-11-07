//Copyright (c) 2014 - 2021, The Trustees of Indiana University.
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

#ifndef GUI_UTIL_COMMAND_HPP
#define GUI_UTIL_COMMAND_HPP

#include <vector>
#include <map>
#include <string>

#include "topfd/common/topfd_para.hpp"

namespace toppic {

namespace command {

std::string geneTopfdCommand(TopfdParaPtr para_ptr,
                             std::vector<std::string> spec_file_lst);  

std::string geneTopIndexCommand(std::map<std::string, std::string> arguments_); 

std::string geneToppicCommand(std::map<std::string, std::string> arguments_, 
                              std::vector<std::string> spec_file_lst_); 

std::string geneTopmgCommand(std::map<std::string, std::string> arguments_, 
                              std::vector<std::string> spec_file_lst_); 

std::string geneTopDiffCommand(std::map<std::string, std::string> arguments_, 
                               std::vector<std::string> spec_file_lst_); 

}
}
#endif
