/*
 * author  Xiaowen Liu
 * date    2013-11-01
 */

#ifndef PROT_COMPONENT_HPP_
#define PROT_COMPONENT_HPP_

#include "acid.hpp"
#include "ptm.hpp"
#include "residue.hpp"

namespace prot {

class Component {
 public:
  Component (const char* acid_file_name,
             const char* ptm_file_name,
             const char* residue_file_name);
  const AcidPtrVec& getAcidList() {return acid_list_;}
  const PtmPtrVec& getPtmList() {return ptm_list_;}
  const ResiduePtrVec& getResidueList() {return residue_list_;}

 private:
  AcidPtrVec acid_list_;
  PtmPtrVec ptm_list_;
  ResiduePtrVec residue_list_;
};

}
#endif

