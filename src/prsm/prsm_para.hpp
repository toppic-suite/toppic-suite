#ifndef PROT_PRSM_PARA_HPP_
#define PROT_PRSM_PARA_HPP_

#include <string>
#include <map>

#include "base/residue.hpp"
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
  const ResiduePtrVec& getFixModResiduePtrVec() {return fix_mod_residue_list_;}
  const ProtModPtrVec& getAllowProtModPtrVec() {return allow_prot_mod_list_;}
  SpParaPtr getSpParaPtr() {return sp_para_ptr_;}

 private:
  std::string search_db_file_name_;
  std::string spec_file_name_;
  std::string log_file_name_;
  std::string exe_dir_;
  int errorTolerance_;

  ResiduePtrVec fix_mod_residue_list_;
  ProtModPtrVec allow_prot_mod_list_;

  /** spectrum parameters */
  SpParaPtr sp_para_ptr_;
};

typedef std::shared_ptr<PrsmPara> PrsmParaPtr;

} /* namespace prot */

#endif /* PRSM_PARA_HPP_ */
