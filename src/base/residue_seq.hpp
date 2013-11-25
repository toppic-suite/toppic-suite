#ifndef PROT_RESIDUE_SEQ_H_
#define PROT_RESIDUE_SEQ_H_

#include "mass_constant.hpp"
#include "residue.hpp"

namespace prot {

class ResidueSeq {
 public:
  ResidueSeq(std::string name, ResiduePtrVec residues);

  /**
   * Returns a sub-peptide of the original peptide.
   **/
  ResidueSeq getSubResidueSeq(int bgn, int end);

  /** Gets length */
  int getLen() {return residues_.size();}

  /** Gets sequences */
  std::string getName() {return name_;}

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

 private:
  /** sequence name */
  std::string name_;
  /** residue list */
  ResiduePtrVec residues_;
  /** the sum of residue mass */
  double residue_mass_sum_;
};

ResidueSeq getEmptyResidueSeq();

}
#endif
