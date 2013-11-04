#include "acid_vec.hpp"

namespace proteomics {
    
	Acid** AcidVec::convert(std::vector<Acid> &acid_vec, const std::string &seq) {
		if (seq.length() == 0) {
            Acid** acids = new Acid*[0];
            return acids;
		} else {
			Acid** acids = new Acid*[seq.length()];
			for (int i = 0; i < seq.length(); i++) {
				acids[i] = AcidVec::getAcidByOneLetter(acid_vec, seq.substr(i, 1));
			}
			return acids;
		}
	}
	
	/**
	 * Returns an amino acid based on the the name. Returns null if the amino
	 * acid name does not exist.
	 */
    Acid* AcidVec::getAcidByName(std::vector<Acid> &acid_vec, const std::string &name) {
		for (int i = 0; i < acid_vec.size(); i++) {
            std::string n = acid_vec[i].getName();
			if (n == name) {
				return &acid_vec[i];
			}
		}
		return NULL;
	}

	/**
	 * Returns an amino acid based on the one letter representation. Returns
	 * null if the one letter representation does not exist.
	 */
	Acid* AcidVec::getAcidByOneLetter(std::vector<Acid> &acid_vec, const std::string &one_letter) {
		for (int i = 0; i < acid_vec.size(); i++) {
            std::string l = acid_vec[i].getOneLetter();
			if (l == one_letter) {
				return &acid_vec[i];
			}
		}
		//logger.debug("Acid not found " + one_letter);
		return NULL;
	}

	/**
	 * Returns an amino acid based on the three letter representation. Returns
	 * null if the three letter representation does not exist.
	 */
	Acid* AcidVec::getAcidByThreeLetter(std::vector<Acid> &acid_vec, const std::string &three_letter) {
		for (int i = 0; i < acid_vec.size(); i++) {
            std::string l = acid_vec[i].getThreeLetter();
			if (l == three_letter) {
				return &acid_vec[i];
			}
		}
		//logger.debug("Acid not found " + three_letter);
		return NULL;
	}
	
};


