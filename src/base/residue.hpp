/*
 * author  Xiaowen Liu
 * date    2013-11-01
 */

#ifndef PROT_BASE_RESIDUE_HPP_
#define PROT_BASE_RESIDUE_HPP_

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

  Residue(xercesc::DOMElement* element); 

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
    return acid_ptr_ == acid_ptr && ptm_ptr_ == ptm_ptr;
  }

  /** Get string representation */
  std::string toString(const std::string &delim_bgn, 
                       const std::string &delim_end);

  std::string toString() {return toString("[", "]");}

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  static std::string getXmlElementName() {return "residue";}

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

}

#endif
