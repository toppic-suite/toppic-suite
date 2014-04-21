/*
 * author  Xiaowen Liu
 * date    2013-11-01
 */

#ifndef PROT_BASE_DATA_HPP_
#define PROT_BASE_DATA_HPP_

#include <map>

#include "base/residue.hpp"
#include "base/prot_mod.hpp"
#include "base/ion_type.hpp"
#include "base/neutral_loss.hpp"
#include "base/activation.hpp"
#include "base/support_peak_type.hpp"

namespace prot {

class BaseData {
 public:
  BaseData (std::string config_file_name);
  BaseData (std::string config_file_name,std::map<std::string,std::string> arguments);

  ResiduePtrVec getFixModResiduePtrVec() {return fix_mod_residue_list_;}
  ProtModPtrVec getAllowProtModPtrVec() {return allow_prot_mod_list_;}
  ActivationPtr getActivationPtr() {return activation_ptr_;}

 private:
  ResiduePtrVec fix_mod_residue_list_;
  ProtModPtrVec allow_prot_mod_list_;
  /* if activation ptr is null, activation types in file are used */
  ActivationPtr activation_ptr_;
};

typedef std::shared_ptr<BaseData> BaseDataPtr;

}
#endif

