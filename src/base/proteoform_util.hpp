#ifndef PROT_BASE_PROTEOFORM_UTIL_HPP_
#define PROT_BASE_PROTEOFORM_UTIL_HPP_

#include "base/proteoform.hpp"


namespace prot {

class ProteoformUtil {
  /* calculate frequencies for n_terminal_residues */
  static ResFreqPtrVec compNTermResidueFreq(const ProteoformPtrVec &prot_mod_forms);

  /* calculater frequences for all residues */
  static ResFreqPtrVec compResidueFreq(const ResiduePtrVec &residue_list,
                                       const ProteoformPtrVec &raw_mods);

  static bool isSameSeqAndMass(ProteoformPtr a, ProteoformPtr b, double ppo);

  static bool isStrictCompatiablePtmSpecies(ProteoformPtr a, ProteoformPtr b, double ppo);

  static ProteoformPtrVec2D divideProteoIntoBlocks(const ProteoformPtrVec &proteo_ptrs, 
                                                   int db_block_size);

};

} /* namespace prot */

#endif 
