#include "residue_seq.hpp"

namespace prot {

ResidueSeq::ResidueSeq(std::string name, ResiduePtrVec residues) {
  name_ = name;
  residues_ = residues_;
  /* get residue mass sum */
  residue_mass_sum_ = 0;
  for (unsigned int i = 0; i < residues_.size(); i++) {
    residue_mass_sum_ += residues_[i]->getMass();
  }
}
ResidueSeq::getSubResidueSeq(int bgn, int end) {
  if (end - bgn < 0) {
    return getEmptyResidueSeq();
  } else {
    ResiduePtrVec sub_residues; 
    std::copy ( v1.begin() + 4, v1.begin() + 8, std::back_inserter(subvector) );
    for (int i = bgn; i <= end; i++) {
      sub_residues.push_back(residues_[i]);
    }
    return ResidueSeq("", subResidues);
  }
}

}
