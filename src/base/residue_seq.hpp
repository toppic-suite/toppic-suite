#ifndef PROT_BASE_RESIDUE_SEQ_HPP_
#define PROT_BASE_RESIDUE_SEQ_HPP_

#include <string>

#include "base/mass_constant.hpp"
#include "base/residue.hpp"

namespace prot {

class ResidueSeq;

typedef std::shared_ptr<ResidueSeq> ResSeqPtr;
typedef std::vector<ResSeqPtr> ResSeqPtrVec;

class ResidueSeq {
 public:
  ResidueSeq(const ResiduePtrVec &residues);

  /**
   * Returns a sub-peptide of the original peptide.
   **/
  ResSeqPtr getSubResidueSeq(int bgn, int end);

  /** Gets length */
  int getLen() {return residues_.size();}

  /** Gets residue at position i */
  ResiduePtr getResiduePtr(int i) {return residues_[i];}

  /** Gets all residues */
  const ResiduePtrVec& getResidues() {return residues_;}

  /** Gets sequence molecular mass */
  double getSeqMass() {
    return residue_mass_sum_ + MassConstant::getWaterMass();
  }

  /** Gets the sum of residue masses */
  double getResMassSum() {return residue_mass_sum_;}

  std::string toString();

  std::string toAcidString();

  static std::string getXmlElementName() {return "residue_seq";}

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

 private:
  /** residue list */
  ResiduePtrVec residues_;
  /** the sum of residue mass */
  double residue_mass_sum_;

  static ResSeqPtr getEmptyResidueSeq();
};


}
#endif
