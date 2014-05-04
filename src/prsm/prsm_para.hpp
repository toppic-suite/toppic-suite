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
  PrsmPara(std::map<std::string,std::string> arguments);
  std::string getSearchDbFileName() {return search_db_file_name_;}
  std::string getSpectrumFileName() {return spectrum_file_name_;}
  ResiduePtrVec getFixModResiduePtrVec() {return fix_mod_residue_list_;}
  ProtModPtrVec getAllowProtModPtrVec() {return allow_prot_mod_list_;}
  SpParaPtr getSpParaPtr() {return sp_para_ptr_;}

 private:
  std::string search_db_file_name_;
  std::string spectrum_file_name_;

  ResiduePtrVec fix_mod_residue_list_;
  ProtModPtrVec allow_prot_mod_list_;

  /** spectrum parameters */
  SpParaPtr sp_para_ptr_;
};

typedef std::shared_ptr<PrsmPara> PrsmParaPtr;

} /* namespace prot */

#endif /* PRSM_PARA_HPP_ */
