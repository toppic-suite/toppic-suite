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

#define SEMI_ALIGN_TYPE_COMPLETE 0
#define SEMI_ALIGN_TYPE_PREFIX   1
#define SEMI_ALIGN_TYPE_SUFFIX   2
#define SEMI_ALIGN_TYPE_INTERNAL 3

#define NEUTRAL_LOSS_NONE "NONE"

class BaseData {
 public:
  BaseData (std::string config_file_name);

  AcidPtrVec& getAcidPtrVec() {return acid_list_;}
  PtmPtrVec& getPtmPtrVec() {return ptm_list_;}
  ResiduePtrVec& getResiduePtrVec() {return residue_list_;}
  TruncPtrVec& getTruncPtrVec() {return trunc_list_;}
  ProtModPtrVec& getProtModPtrVec() {return prot_mod_list_;}
  IonTypePtrVec& getIonTypePtrVec() {return ion_type_list_;}
  NeutralLossPtrVec& getNeutralLossPtrVec() {return neutral_loss_list_;}
  NeutralLossPtr getNeutralLossNonePtr() {
    return getNeutralLossPtrByName(neutral_loss_list_, NEUTRAL_LOSS_NONE);
  }
  ActivationPtrVec& getActivationPtrVec() {return activation_list_;}

  ResiduePtrVec& getFixModResiduePtrVec() {return fix_mod_residue_list_;}
  ProtModPtrVec& getAllowProtModPtrVec() {return allow_prot_mod_list_;}
  ActivationPtr& getActivationPtr() {return activation_ptr_;}

 private:
  AcidPtrVec acid_list_;
  PtmPtrVec ptm_list_;
  ResiduePtrVec residue_list_;

  TruncPtrVec trunc_list_;
  ProtModPtrVec prot_mod_list_;

  IonTypePtrVec ion_type_list_;
  NeutralLossPtrVec neutral_loss_list_;

  ActivationPtrVec activation_list_;

  /* configuration */
  ResiduePtrVec fix_mod_residue_list_;
  ProtModPtrVec allow_prot_mod_list_;
  /* if activation ptr is null, activation types in file are used */
  ActivationPtr activation_ptr_;
};

typedef std::shared_ptr<BaseData> BaseDataPtr;

}
#endif

