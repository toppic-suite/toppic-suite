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

#ifndef TOPPIC_STAT_MCMC_MASS_TABLE_UTIL_HPP_
#define TOPPIC_STAT_MCMC_MASS_TABLE_UTIL_HPP_

#include <vector>
#include <set>
#include <map>
#include <string>

#include "common/base/residue.hpp"

#include "para/sp_para.hpp"

#include "stat/mcmc/mcmc_mng.hpp"

namespace toppic {
namespace mass_table_util {

std::map<int, std::vector<std::string> > readMassTable(MCMCMngPtr mng_ptr);

std::map<int, std::vector<std::string> > geneMassTable(MCMCMngPtr mng_ptr);

std::map<int, std::vector<std::string> > geneMassTableFixMod(MCMCMngPtr mng_ptr);

}  // namespace mass_table_util
}  // namespace toppic

#endif
