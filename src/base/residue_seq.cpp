#include "residue_seq.hpp"

namespace prot {

ResidueSeq::ResidueSeq(ResiduePtrVec residues) {
  residues_ = residues;
  /* get residue mass sum */
  residue_mass_sum_ = 0;
  for (unsigned int i = 0; i < residues_.size(); i++) {
    residue_mass_sum_ += residues_[i]->getMass();
  }
}


ResidueSeq ResidueSeq::getSubResidueSeq(int bgn, int end) {
  if (end - bgn < 0) {
    return getEmptyResidueSeq();
  } else {
    ResiduePtrVec sub_residues; 
    std::copy (residues_.begin() + bgn, residues_.begin() + end, 
               std::back_inserter(sub_residues) );
    return ResidueSeq(sub_residues);
  }
}

ResidueSeq getEmptyResidueSeq() {
  ResiduePtrVec residues;
  return ResidueSeq(residues);
}

}
