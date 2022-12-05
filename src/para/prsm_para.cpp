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

#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"
#include "common/util/file_util.hpp"

#include "common/base/mod_util.hpp"
#include "common/base/prot_mod_base.hpp"
#include "para/prsm_para.hpp"

namespace toppic {

PrsmPara::PrsmPara(std::map<std::string, std::string> &arguments) {
  ori_db_name_ = arguments["oriDatabaseFileName"];

  search_db_file_name_ = file_util::filenameFromEntirePath(arguments["databaseFileName"]);

  search_db_file_name_with_folder_ = ori_db_name_ + "_idx" 
    + file_util::getFileSeparator() + search_db_file_name_;

  spec_file_name_ = arguments["spectrumFileName"];

  resource_dir_ = arguments["resourceDir"];

  group_spec_num_ = std::stoi(arguments["groupSpectrumNumber"]);

  fix_mod_list_ = mod_util::geneFixedModList(arguments["fixedMod"]);

  std::string prot_mod_str = arguments["allowProtMod"];
  std::vector<std::string> strs = str_util::split(prot_mod_str, ",");
  for (size_t i = 0; i < strs.size(); i++) {
    ProtModPtrVec mods = ProtModBase::getProtModPtrByType(strs[i]);
    LOG_DEBUG("prot mod type " << strs[i] << " num " << mods.size());
    prot_mod_list_.insert(prot_mod_list_.end(), mods.begin(), mods.end());
  }

  std::string activation_name = arguments["activation"];
  int ppm = std::stoi(arguments["massErrorTolerance"]);
  double n_term_label_mass = 0;
  if (arguments["nTermLabelMass"] != "") {
    n_term_label_mass = std::stod(arguments["nTermLabelMass"]);
  }
  sp_para_ptr_ = std::make_shared<SpPara>(activation_name, n_term_label_mass, ppm); 
}

} /* namespace toppic */
