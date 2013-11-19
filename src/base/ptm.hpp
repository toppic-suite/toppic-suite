/*
 * author  Xiaowen Liu
 * date    2013-11-17
 */

#ifndef PROTOMICS_PTM_H_
#define PROTOMICS_PTM_H_

#include <string>

namespace proteomics {
class Ptm {
    public:
    Ptm(const std::string &abbr_name, 
            const Acid *valid_acids[], 
            double mono_mass, 
            bool is_empty); 
    ~Ptm();
    
	/* Get amino acid composition. */
    std::string getAbbrName() { return abbr_name_;}

	/* Get  monoisotopic mass. */
	double getMonoMass() {return mono_mass_;}

	/* Get valid acid list. */
    Acid* getValidAcids() {return valid_acids_;}

	/* Is it an empty PTM. */
    bool isEmpty() {return is_empty_;}

    private:
	/* Abbreviation name of a PTM */
    std::string abbr_name_;
	/* Valid acids of the PTM */
    Acid valid_acids_[];
	/* Is it an empty PTM */
    bool is_empty_;
};
}
#endif

