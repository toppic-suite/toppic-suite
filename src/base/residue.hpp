/*
 * author  Xiaowen Liu
 * date    2013-11-01
 */

#ifndef PROT_RESIDUE_HPP_
#define PROT_RESIDUE_HPP_

#include <string>
#include <memory>
#include <map>

#include "base/acid.hpp"
#include "base/ptm.hpp"
#include "base/logger.hpp"

namespace prot {

class Residue {
 public:
  Residue(AcidPtr acid_ptr, PtmPtr ptm_ptr); 

  Residue(const std::string &acid_name, const std::string &abbr_name);

  /** Get amino acid. */
  AcidPtr getAcidPtr() {return acid_ptr_; }

  /** Get residue mass. */
  double getMass() { return mass_; }

  /** Get post-translational modification. */
  PtmPtr getPtmPtr() { return ptm_ptr_; }

  /**
   * Checks if the residue contains the same amino acid and ptm.
   */
  bool isSame(AcidPtr acid_ptr, PtmPtr ptm_ptr) {
    return acid_ptr_.get() == acid_ptr.get() && ptm_ptr_.get() == ptm_ptr.get();
  }

  /** Get string representation */
  std::string toString(const std::string &delim_bgn, 
                       const std::string &delim_end);

  std::string toString() {return toString("[", "]");}

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

private:
  /** amino acid */
  AcidPtr acid_ptr_;
  /** post-translational modification */
  PtmPtr ptm_ptr_;
  /** residue mass */
  double mass_;
};

typedef std::shared_ptr<Residue> ResiduePtr;
typedef std::vector<ResiduePtr> ResiduePtrVec;
typedef std::vector<ResiduePtrVec> ResiduePtrVec2D;

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
