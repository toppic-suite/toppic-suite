/*
 * author  Xiaowen Liu
 * date    2013-11-01
 */

#ifndef PROTOMICS_COMPONENT_H_
#define PROTOMICS_COMPONENT_H_

#include "acid.hpp"
#include "ptm.hpp"
#include "residue.hpp"

namespace proteomics {

class Component {
 public:
  Component (const char* acid_file_name,
             const char* ptm_file_name,
             const char* residue_file_name);
  const AcidPtrVec& getAcidPtrVec() {return acid_ptr_vec_;}
  const PtmPtrVec& getPtrPtrVec() {return ptm_ptr_vec_;}
  const ResiduePtrVec& getResiduePtrVec() {return residue_ptr_vec_;}

 private:
  AcidPtrVec acid_ptr_vec_;
  PtmPtrVec ptm_ptr_vec_;
  ResiduePtrVec residue_ptr_vec_;
};

}
#endif

