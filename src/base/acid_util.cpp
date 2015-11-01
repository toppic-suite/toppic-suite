#include "base/logger.hpp"
#include "base/acid_util.hpp"

namespace prot {

/**
 * Converts a protein sequence (with one letter representation of amino
 * acids) to an amino acid array.
 */
AcidPtrVec AcidUtil::convertStrToAcidPtrVec(const std::string &seq) {
  AcidPtrVec acid_seq;
  for (size_t i = 0; i < seq.length(); i++) {
    acid_seq.push_back(AcidBase::getAcidPtrByOneLetter(seq.substr(i, 1)));
  }
  return acid_seq;
}

double AcidUtil::compAcidPtrVecMass(const AcidPtrVec &ptr_vec) {
  double mass = 0;
  for (size_t i = 0; i < ptr_vec.size(); i++) {
    mass += ptr_vec[i]->getMonoMass();
  }
  return mass;
}

} /* end namespace */

