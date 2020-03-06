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

#include <fstream>

#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "common/base/activation_base.hpp"
#include "common/base/mass_constant.hpp"
#include "common/base/mod_util.hpp"
#include "common/base/prot_mod_base.hpp"
#include "para/prsm_para.hpp"

namespace toppic {

PrsmPara::PrsmPara(std::map<std::string, std::string> &arguments) {
  search_db_file_name_ = arguments["databaseFileName"];

  ori_db_name_ = arguments["oriDatabaseFileName"];

  spec_file_name_ = arguments["spectrumFileName"];

  resource_dir_ = arguments["resourceDir"];

  ppm_ = std::stoi(arguments["massErrorTolerance"]);

  group_spec_num_ = std::stoi(arguments["groupSpectrumNumber"]);

  fix_mod_list_ = mod_util::geneFixedModList(arguments["fixedMod"]);

  activation_ = arguments["activation"];

  fixed_mod_ = arguments["fixedMod"];

  error_tol_ = arguments["massErrorTolerance"];

  allow_prot_mod_ = arguments["allowProtMod"];

  decoy_ = arguments["searchType"];

  std::string prot_mod_str = arguments["allowProtMod"];
  //boost::split(strs, prot_mod_str, boost::is_any_of(","));
  std::vector<std::string> strs = str_util::split(prot_mod_str, ",");
  for (size_t i = 0; i < strs.size(); i++) {
    ProtModPtrVec mods = ProtModBase::getProtModPtrByType(strs[i]);
    LOG_DEBUG("prot mod type " << strs[i] << " num " << mods.size());
    prot_mod_list_.insert(prot_mod_list_.end(), mods.begin(), mods.end());
  }

  std::string activation_name = arguments["activation"];

  std::set<std::string> skip_list;

  if (arguments["skipList"] != "") {
    std::ifstream infile(arguments["skipList"]);
    std::string line;
    while (std::getline(infile, line)) {
      //boost::split(strs, line, boost::is_any_of(" "));
      std::vector<std::string> strs = str_util::split(line, " "); 
      for (size_t i = 0; i < strs.size(); i++) {
        skip_list.insert(strs[i]);
      }
    }
    infile.close();
  }

  sp_para_ptr_ = std::make_shared<SpPara>(activation_name, ppm_); 
  
  sp_para_ptr_->setSkipList(skip_list);
}

} /* namespace toppic */
