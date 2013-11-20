#ifndef PROT_ACID_UTIL_H_
#define PROT_ACID_UTIL_H_

#include <vector>

#include "acid.hpp"

namespace prot {


AcidPtrVec getAcidPtrVecInstance(const char* file_name);

/**
 * Returns an amino acid based on the the name. Returns null if the amino
 * acid name does not exist.
 */
AcidPtr getAcidPtrByName(AcidPtrVec &acid_ptr_vec, 
                         const std::string &name);
/**
 *
 * Returns an amino acid based on the one letter representation. Returns
 * null if the one letter representation does not exist.
 */
AcidPtr getAcidPtrByOneLetter(AcidPtrVec &acid_ptr_vec, 
                              const std::string &one_letter);
/**
 * Returns an amino acid based on the three letter representation. Returns
 * null if the three letter representation does not exist.
 */
AcidPtr getAcidPtrByThreeLetter(AcidPtrVec &acid_ptr_vec, 
                                const std::string &three_letter);

/**
 * Checks if the list contains an amino acid with the specific name.
 */
bool containsName(AcidPtrVec &acid_ptr_vec, const std::string &name);

/**
 * Checks if the list contains an amino acid with the specific one letter
 * representation.
 */
bool containsOneLetter(AcidPtrVec &acid_ptr_vec, const std::string &one_letter);

/**
 * Checks if the list contains an amino acid with the specific three letter
 * representation.
 */
bool containsThreeLetter(AcidPtrVec &acid_ptr_vec, const std::string &three_letter);

/**
 * Converts a protein sequence (with one letter representation of amino
 * acids) to an amino acid array.
 */
AcidPtrVec convertSeqToAcidPtrVec(AcidPtrVec &acid_ptr_vec, 
                                  const std::string &seq);
}
#endif
