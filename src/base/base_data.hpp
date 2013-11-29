/*
 * author  Xiaowen Liu
 * date    2013-11-01
 */

#ifndef PROT_BASE_DATA_HPP_
#define PROT_BASE_DATA_HPP_

#include "acid.hpp"
#include "ptm.hpp"
#include "residue.hpp"

namespace prot {

class BaseData {
 public:
  BaseData (const char* config_file_name);

  const AcidPtrVec& getAcidPtrVec() {return acid_list_;}
  const PtmPtrVec& getPtmPtrVec() {return ptm_list_;}
  const ResiduePtrVec& getResiduePtrVec() {return residue_list_;}

 private:
  AcidPtrVec acid_list_;
  PtmPtrVec ptm_list_;
  ResiduePtrVec residue_list_;
};

}
#endif

