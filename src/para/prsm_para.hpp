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

#ifndef TOPPIC_PARA_PRSM_PARA_HPP_
#define TOPPIC_PARA_PRSM_PARA_HPP_

#include <map>

#include "common/base/mod.hpp"
#include "common/base/prot_mod.hpp"
#include "para/sp_para.hpp"

namespace toppic {

class PrsmPara {
 public:
  PrsmPara(std::map<std::string,std::string> &arguments);

  std::string getSearchDbFileName() {return search_db_file_name_;}

  std::string getOriDbName(){return ori_db_name_;}
  
  std::string getSpectrumFileName() {return spec_file_name_;}

  std::string getResourceDir() {return resource_dir_;}

  std::string getActivation() {return activation_;}

  std::string getFixedMod() {return fixed_mod_;}

  std::string getErrorTolerance() {return error_tol_;}

  std::string getProtMod() {return allow_prot_mod_;}

  std::string getSearchType() {return decoy_;}

  int getPPM() {return ppm_;}

  int getGroupSpecNum() {return group_spec_num_;}

  const ModPtrVec& getFixModPtrVec() {return fix_mod_list_;}

  const ProtModPtrVec& getProtModPtrVec() {return prot_mod_list_;}

  SpParaPtr getSpParaPtr() {return sp_para_ptr_;}

  bool doLocaliztion() {return localization_;}

  std::string getIndexDir() {return index_file_dir_;}

  void setIndexDir(std::string dir) {index_file_dir_ = dir;}

 private:
  std::string search_db_file_name_;

  std::string ori_db_name_;

  std::string spec_file_name_;

  std::string resource_dir_;

  std::string activation_;

  std::string fixed_mod_;

  std::string error_tol_;

  std::string allow_prot_mod_;

  std::string decoy_;

  int ppm_;

  ModPtrVec fix_mod_list_;

  ProtModPtrVec prot_mod_list_;

  int group_spec_num_;

  bool localization_;

  std::string index_file_dir_;

  /** spectrum parameters */
  SpParaPtr sp_para_ptr_;
};

typedef std::shared_ptr<PrsmPara> PrsmParaPtr;

} /* namespace toppic */

#endif /* PRSM_PARA_HPP_ */
