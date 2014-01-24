/*
 * author  Xiaowen Liu
 * date    2013-11-01
 */

#ifndef PROT_BASE_DATA_HPP_
#define PROT_BASE_DATA_HPP_

#include "base/acid.hpp"
#include "base/ptm.hpp"
#include "base/residue.hpp"
#include "base/trunc.hpp"
#include "base/prot_mod.hpp"
#include "base/ion_type.hpp"
#include "base/neutral_loss.hpp"
#include "base/activation.hpp"

namespace prot {

#define DEFAULT_TYPE_NAME "NONE"
#define ACETYLATION_NAME "ACETYLATION"

class BaseData {
 public:
  BaseData (std::string config_file_name);

  //TruncPtrVec& getTruncPtrVec() {return trunc_list_;}
  ProtModPtrVec& getProtModPtrVec() {return prot_mod_list_;}
  ProtModPtr getDefaultProtModPtr() {
    return getProtModPtrByName(prot_mod_list_, DEFAULT_TYPE_NAME);
  }
  ProtModPtr getAcetylationProtModPtr() {
    return getProtModPtrByName(prot_mod_list_, ACETYLATION_NAME);
  }
  //IonTypePtrVec& getIonTypePtrVec() {return ion_type_list_;}
  NeutralLossPtrVec& getNeutralLossPtrVec() {return neutral_loss_list_;}
  NeutralLossPtr getDefaultNeutralLossPtr() {
    return getNeutralLossPtrByName(neutral_loss_list_, DEFAULT_TYPE_NAME);
  }
  //ActivationPtrVec& getActivationPtrVec() {return activation_list_;}

  ResiduePtrVec& getFixModResiduePtrVec() {return fix_mod_residue_list_;}
  ProtModPtrVec& getAllowProtModPtrVec() {return allow_prot_mod_list_;}
  ActivationPtr& getActivationPtr() {return activation_ptr_;}

 private:
  //TruncPtrVec trunc_list_;
  ProtModPtrVec prot_mod_list_;

  //IonTypePtrVec ion_type_list_;
  NeutralLossPtrVec neutral_loss_list_;

  //ActivationPtrVec activation_list_;

  /* configuration */
  ResiduePtrVec fix_mod_residue_list_;
  ProtModPtrVec allow_prot_mod_list_;
  /* if activation ptr is null, activation types in file are used */
  ActivationPtr activation_ptr_;
};

typedef std::shared_ptr<BaseData> BaseDataPtr;

}
#endif

