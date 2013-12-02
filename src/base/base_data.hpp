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

namespace prot {

class BaseData {
 public:
  BaseData (char const* config_file_name);

  const AcidPtrVec& getAcidPtrVec() {return acid_list_;}
  PtmPtrVec& getPtmPtrVec() {return ptm_list_;}
  ResiduePtrVec& getResiduePtrVec() {return residue_list_;}
  TruncPtrVec& getTruncPtrVec() {return trunc_list_;}
  ProtModPtrVec& getProtModPtrVec() {return prot_mod_list_;}
  IonTypePtrVec& getIonTypePtrVec() {return ion_type_list_;}
  NeutralLossPtrVec& getNeutralLossPtrVec() {return neutral_loss_list_;}

 private:
  AcidPtrVec acid_list_;
  PtmPtrVec ptm_list_;
  ResiduePtrVec residue_list_;

  TruncPtrVec trunc_list_;
  ProtModPtrVec prot_mod_list_;

  IonTypePtrVec ion_type_list_;
  NeutralLossPtrVec neutral_loss_list_;
};

typedef std::shared_ptr<BaseData> BaseDataPtr;

}
#endif

