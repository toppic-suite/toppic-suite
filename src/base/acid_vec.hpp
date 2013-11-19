#ifndef PROTOMICS_ACID_VEC_H_
#define PROTOMICS_ACID_VEC_H_

#include <vector>

#include "acid.hpp"

namespace proteomics {
class AcidVec {
public:
	static std::vector<Acid> getInstance();
    
	/**
	 * Returns an amino acid based on the the name. Returns null if the amino
	 * acid name does not exist.
	 */
	static Acid getAcidByName(std::vector<Acid> &acid_vec, const std::string &name);
	/**
     *
	 * Returns an amino acid based on the one letter representation. Returns
	 * null if the one letter representation does not exist.
	 */
	static Acid getAcidByOneLetter(std::vector<Acid> &acid_vec, const std::string &one_letter);
	/**
	 * Returns an amino acid based on the three letter representation. Returns
	 * null if the three letter representation does not exist.
	 */
	static Acid getAcidByThreeLetter(std::vector<Acid> &acid_vec, const std::string &three_letter);

	/**
	 * Checks if the list contains an amino acid with the specific name.
	 */
	static bool containsName(std::vector<Acid> &acid_vec, const std::string &name);

	/**
	 * Checks if the list contains an amino acid with the specific one letter
	 * representation.
	 */
	static bool containsOneLetter(std::vector<Acid> &acid_vec, const std::string &one_letter);

	/**
	 * Checks if the list contains an amino acid with the specific three letter
	 * representation.
	 */
	static bool containsThreeLetter(std::vector<Acid> &acid_vec, const std::string &three_letter);

	/**
	 * Converts a protein sequence (with one letter representation of amino
	 * acids) to an amino acid array.
	 */
	static std::vector<Acid> convert(std::vector<Acid> &acid_vec, const std::string &seq);
};
}
#endif
