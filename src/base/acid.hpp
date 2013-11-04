/*
 * author  Xiaowen Liu
 * date    2013-11-01
 */

#ifndef PROTOMICS_ACID_H_
#define PROTOMICS_ACID_H_

#include <string>

#include <xercesc/dom/DOM.hpp>

namespace proteomics {
class Acid {
    public:
    Acid (const std::string &name, const std::string &one_letter, const std::string &three_letter,
            const std::string &composition, double mono_mass, double avg_mass); 
    
    Acid (xercesc::DOMElement * element);

	/* Get amino acid composition. */
    std::string get_acid_composition() { return composition_;}

	/* Get average mass. */
	double get_avg_mass() {return avg_mass_;}

	/* Get  monoisotopic mass. */
	double get_mono_mass() {return mono_mass_;}

	/* Get amino acid name. */
    std::string get_name() {return name_;}

	/* Get amino acid one letter representation. */
    std::string get_one_letter() {return one_letter_;}

	/* Get amino acid three letter representation. */
    std::string get_three_letter() {return three_letter_;}

    private:
	/* Name of amino acid */
    std::string name_;
	/* One letter representation */
    std::string one_letter_;
	/* Three letter representation */
    std::string three_letter_;
	/* amino acid chemical composition */
    std::string composition_;
	/* residue monoisotopic mass */
    double mono_mass_;
	/* residue average mass */
    double avg_mass_;
};
}
#endif

