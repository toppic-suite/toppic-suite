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


#include "spec/ms_header.hpp"
#include "deconv/env/match_env.hpp"

namespace toppic {

namespace match_env_writer {

void write_env(std::ofstream &file, MsHeaderPtr header, MatchEnvPtr match_env); 

void write_env_vec(std::ofstream &file, MsHeaderPtr header, const MatchEnvPtrVec & envs); 

void write(const std::string & file, MsHeaderPtr header, const MatchEnvPtrVec & envs); 


}  // namespace msalign_writer

}  // namespace toppic
