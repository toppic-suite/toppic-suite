//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
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

#ifndef TOPPIC_CONSOLE_TOPFD_PROCESS_HPP_
#define TOPPIC_CONSOLE_TOPFD_PROCESS_HPP_

#include <string>
#include <map>

namespace toppic {

namespace topfd_process {

std::string geneArgumentStr(std::map<std::string, std::string> arguments, 
                            const std::string & prefix);

int processOneFile(std::map<std::string, std::string> arguments, 
                   const std::string &argument_str,
                   const std::string &spec_file_name, int frac_id);

int process(std::map<std::string, std::string> arguments, 
            std::vector<std::string> spec_file_lst); 

}

}

#endif
