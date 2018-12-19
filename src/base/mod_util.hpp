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


#ifndef TOPPIC_BASE_MOD_UTIL_HPP_
#define TOPPIC_BASE_MOD_UTIL_HPP_

#include <string>
#include <vector>

#include "base/mod.hpp"
#include "base/residue.hpp"

namespace toppic {

namespace mod_util {

ModPtrVec readModXml(const std::string &file_name);

std::vector<ModPtrVec> readModTxt(const std::string &file_name);

ModPtrVec geneFixedModList(const std::string &str);

ResiduePtrVec geneResidueListWithMod(const ResiduePtrVec & residue_list,
                                     const ModPtrVec & fix_mod_list);

std::vector<double> getModMassVec(const ModPtrVec & var_mod_list);

}  // namespace mod_util

}  // namespace toppic

#endif
