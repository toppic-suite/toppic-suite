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


//#include <iomanip>

#include "spec/ms_header.hpp"
#include "feature/match_env.hpp"

namespace prot {

namespace match_env_writer {

void write_env(std::ofstream &file, MatchEnvPtr match_env); 

void write_spectrum(std::ofstream &file, MsHeaderPtr header, 
                    MatchEnvPtrVec &envs); 

void write(MsHeaderPtr header, MatchEnvPtrVec &envs); 

}  // namespace msalign_writer

}  // namespace prot
