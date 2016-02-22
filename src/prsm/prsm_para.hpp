#ifndef PROT_PRSM_PRSM_PARA_HPP_
#define PROT_PRSM_PRSM_PARA_HPP_

#include <string>
#include <map>

#include "base/mod.hpp"
#include "base/prot_mod.hpp"
#include "base/activation.hpp"
#include "spec/peak_tolerance.hpp"
#include "spec/sp_para.hpp"

namespace prot {

class PrsmPara {
 public:
  PrsmPara(std::map<std::string,std::string> &arguments);
  std::string getSearchDbFileName() {return search_db_file_name_;}
  std::string getSpectrumFileName() {return spec_file_name_;}
  std::string getLogFileName() {return log_file_name_;}
  std::string getExeDir(){return exe_dir_;}
  int getErrorTolerance(){return errorTolerance_;}
  int getGroupSpecNum() {return group_spec_num_;}
  const ModPtrVec& getFixModPtrVec() {return fix_mod_list_;}
  const ProtModPtrVec& getProtModPtrVec() {return prot_mod_list_;}
  SpParaPtr getSpParaPtr() {return sp_para_ptr_;}
  bool doLocaliztion() {return localization_;}

 private:
  std::string search_db_file_name_;
  std::string spec_file_name_;
  std::string log_file_name_;
  std::string exe_dir_;
  int errorTolerance_;

  ModPtrVec fix_mod_list_;
  ProtModPtrVec prot_mod_list_;

  int group_spec_num_;

  bool localization_;

  /** spectrum parameters */
  SpParaPtr sp_para_ptr_;
};

typedef std::shared_ptr<PrsmPara> PrsmParaPtr;

} /* namespace prot */

#endif /* PRSM_PARA_HPP_ */
