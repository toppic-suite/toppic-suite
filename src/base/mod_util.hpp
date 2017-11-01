//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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


#ifndef PROT_BASE_MOD_UTIL_HPP_
#define PROT_BASE_MOD_UTIL_HPP_

#include "base/mod.hpp"

namespace prot {

class ModUtil {
 public:
  static ModPtrVec readModXml(const std::string &file_name);

  static std::vector<ModPtrVec> readModTxt(const std::string &file_name);

  static ModPtrVec geneFixedModList(const std::string &str);

  static ResiduePtrVec geneResidueListWithMod(ResiduePtrVec residue_list,
                                              ModPtrVec fix_mod_list);

  static std::vector<double> getModMassVec(const ModPtrVec & var_mod_list);
};

}  // namespace prot
#endif
