#ifndef PROT_RESIDUE_SEQ_HPP_
#define PROT_RESIDUE_SEQ_HPP_

#include <string>

#include "base/mass_constant.hpp"
#include "base/residue.hpp"
#include "base/prot_mod.hpp"

namespace prot {

class ResidueSeq {
 public:
  ResidueSeq(ResiduePtrVec residues);

  /**
   * Returns a sub-peptide of the original peptide.
   **/
  ResidueSeq getSubResidueSeq(int bgn, int end);

  /** Gets length */
  int getLen() {return residues_.size();}

  /** Gets residue at position i */
  ResiduePtr getResiduePtr(int i) {return residues_[i];}

  /** Gets all residues */
  ResiduePtrVec getResidues() {return residues_;}

  /** Gets sequence molecular mass */
  double getSeqMass() {
    return residue_mass_sum_ + MassConstant::getWaterMass();
  }

  /** Gets the sum of residue masses */
  double getResMassSum() {return residue_mass_sum_;}

  std::string toString();

  void appendXml(XmlDOMDocument* xml_doc,xercesc::DOMElement* parent);

  bool allowsMod(ProtModPtr mod);
 private:
  /** residue list */
  ResiduePtrVec residues_;
  /** the sum of residue mass */
  double residue_mass_sum_;
};

typedef std::shared_ptr<ResidueSeq> ResSeqPtr;
typedef std::vector<ResSeqPtr> ResSeqPtrVec;

ResidueSeq getEmptyResidueSeq();

}
#endif
