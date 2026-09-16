// Copyright (c) 2014 - 2026, The Trustees of Indiana University, Tulane
// University.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "para/prsm_para.hpp"

#include <cstddef>
#include <string>
#include <vector>

#include "common/base/mod_util.hpp"
#include "common/base/prot_mod_base.hpp"
#include "common/util/file_util.hpp"
#include "common/util/logger.hpp"
#include "common/util/str_util.hpp"

namespace toppic {

namespace {

// Look up an argument, returning an empty string when the key is absent. This
// preserves the old operator[] behaviour without mutating the (now const) map.
std::string getArg(const std::map<std::string, std::string>& arguments,
                   const std::string& key) {
  auto it = arguments.find(key);
  return it == arguments.end() ? std::string() : it->second;
}

}  // namespace

PrsmPara::PrsmPara(const std::map<std::string, std::string>& arguments) {
  ori_db_name_ = getArg(arguments, "oriDatabaseFileName");

  search_db_file_name_ =
      file_util::filenameFromEntirePath(getArg(arguments, "databaseFileName"));

  search_db_file_name_with_folder_ = ori_db_name_ + "_idx" +
                                     file_util::getFileSeparator() +
                                     search_db_file_name_;

  spec_file_name_ = getArg(arguments, "spectrumFileName");

  resource_dir_ = getArg(arguments, "resourceDir");

  group_spec_num_ = std::stoi(getArg(arguments, "groupSpectrumNumber"));

  fix_mod_list_ = mod_util::geneFixedModList(getArg(arguments, "fixedMod"));

  std::string prot_mod_str = getArg(arguments, "allowProtMod");
  std::vector<std::string> strs = str_util::split(prot_mod_str, ",");
  for (std::size_t i = 0; i < strs.size(); i++) {
    ProtModPtrVec mods = ProtModBase::getProtModPtrByType(strs[i]);
    LOG_DEBUG("prot mod type " << strs[i] << " num " << mods.size());
    prot_mod_list_.insert(prot_mod_list_.end(), mods.begin(), mods.end());
  }

  std::string prot_type_str = getArg(arguments, "allowProtType");
  std::vector<std::string> type_strs = str_util::split(prot_type_str, ",");
  for (std::size_t i = 0; i < type_strs.size(); i++) {
    ProteoformTypePtr prot_type =
        ProteoformType::getProtTypePtrByName(type_strs[i]);
    if (prot_type != nullptr) {
      LOG_DEBUG("Proteoform type " << prot_type->getName());
      prot_type_list_.push_back(prot_type);
    }
  }

  std::string activation_name = getArg(arguments, "activation");
  int ppm = std::stoi(getArg(arguments, "massErrorTolerance"));
  double n_term_label_mass = 0;
  if (getArg(arguments, "nTermLabelMass") != "") {
    n_term_label_mass = std::stod(getArg(arguments, "nTermLabelMass"));
  }
  sp_para_ptr_ =
      std::make_shared<SpPara>(activation_name, n_term_label_mass, ppm);
  if (getArg(arguments, "envCnnCutoff") != "") {
    sp_para_ptr_->setEnvCnnCutoff(std::stod(getArg(arguments, "envCnnCutoff")));
  }
}

bool PrsmPara::allowProtType(const ProteoformTypePtr& type_ptr) const {
  for (std::size_t i = 0; i < prot_type_list_.size(); i++) {
    if (prot_type_list_[i] == type_ptr) {
      return true;
    }
  }

  return false;
}

} /* namespace toppic */
