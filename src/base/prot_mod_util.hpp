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


#ifndef PROT_BASE_PROT_MOD_UTIL_HPP_
#define PROT_BASE_PROT_MOD_UTIL_HPP_

#include <string>

#include "base/prot_mod.hpp"
#include "base/residue.hpp"

namespace prot {

namespace prot_mod_util {

bool allowMod(ProtModPtr prot_mod_ptr, const ResiduePtrVec &residues);

ProtModPtrVec readProtMod(const std::string &file_name);

ProtModPtr findNME_Acetylation(const ProtModPtrVec &prot_mod_ptrs, 
                               const ResiduePtrVec &residues);

}  // namespace prot_mod_util

}  // namespace prot
#endif
