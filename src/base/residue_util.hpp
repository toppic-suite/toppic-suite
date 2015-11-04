/*
 * author  Xiaowen Liu
 * date    2013-11-01
 */

#ifndef PROT_BASE_RESIDUE_BASE_HPP_
#define PROT_BASE_RESIDUE_BASE_HPP_

#include <string>
#include <memory>
#include <map>

#include "base/residue.hpp"
#include "base/logger.hpp"

namespace prot {

ResiduePtr getResiduePtrByAcid(const ResiduePtrVec &residue_list,
                               AcidPtr acid_ptr);

int findResidue(const ResiduePtrVec &residue_list, ResiduePtr residue_ptr);

ResiduePtrVec convertAcidToResidueSeq(const ResiduePtrVec &residue_list,
                                      const AcidPtrVec &acid_list);
/* residue factory */
class ResidueFactory {
 public:
  static void initFactory(const std::string &file_name);

  static const ResiduePtrVec& getBaseResiduePtrVec() {return residue_ptr_vec_;}
  
  static ResiduePtr getBaseResiduePtrByAcidPtm(AcidPtr acid_ptr, PtmPtr ptm_ptr);
  
  static ResiduePtr addBaseResidue(AcidPtr acid_ptr, PtmPtr ptm_ptr);
  
  static ResiduePtrVec getResiduePtrVecInstance(const std::string &file_name);

 private:
  static ResiduePtrVec residue_ptr_vec_;
};

/* residue list factory */
class FixResidueFactory {
 public:
  static void initFactory(const std::string &file_name);

  static ResiduePtrVec getFixResiduePtrVec(const std::string &id);
  
 private:
  static std::map<std::string,ResiduePtrVec> fix_res_list_map_;
};

}
#endif
